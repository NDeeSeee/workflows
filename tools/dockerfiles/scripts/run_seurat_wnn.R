#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB should be good by default

suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(future))
suppressMessages(library(DESeq2))
suppressMessages(library(tibble))
suppressMessages(library(scales))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(argparse))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(reticulate))
suppressMessages(library(sctransform))
suppressMessages(library(rtracklayer))
suppressMessages(library(RColorBrewer))
suppressMessages(library(SeuratWrappers))


set_threads <- function (threads) {
    invisible(capture.output(plan("multiprocess", workers=threads)))
    invisible(capture.output(plan()))
    invisible(capture.output(setDTthreads(threads)))
}


writeSparseMatrix <- function (inMat, outFname, sliceSize=1000){
    fnames <- c()
    mat <- inMat
    geneCount <- nrow(mat)
    startIdx <- 1
    while (startIdx < geneCount){
        endIdx <- min(startIdx + sliceSize - 1, geneCount)
        matSlice <- mat[startIdx:endIdx,]
        denseSlice <- as.matrix(matSlice)
        dt <- data.table(denseSlice)
        dt <- cbind(gene=rownames(matSlice), dt)
        writeHeader <- FALSE
        if (startIdx == 1) { 
            writeHeader <- TRUE
        }
        sliceFname <- paste0("temp", startIdx, ".txt")
        fwrite(dt, sep="\t", file=sliceFname, quote=FALSE, col.names=writeHeader)
        fnames <- append(fnames, sliceFname)
        startIdx <- startIdx + sliceSize
    }
    system(paste("cat", paste(fnames, collapse=" "), "| gzip >", outFname, sep=" "))
    unlink(fnames)
}


findMatrix <- function(object, matrix.slot){
    if (matrix.slot == "counts") {
        counts <- GetAssayData(object=object, slot="counts")
    } else if (matrix.slot == "scale.data") {
        counts <- GetAssayData(object=object, slot="scale.data")
    } else if (matrix.slot == "data") {
        counts <- GetAssayData(object=object)
    } else {
        error("matrix.slot can only be one of: counts, scale.data, data")
    }
}


ExportToCellbrowser <- function(
    object,
    matrix.slot,
    dir,
    cluster.field,
    features,
    meta.fields,
    meta.fields.names
) {
    if (!dir.exists(dir)) {
        dir.create(dir)
    }
    idents <- Idents(object)
    meta <- object@meta.data
    cellOrder <- colnames(object)
    counts <- findMatrix(object, matrix.slot)
    genes <- rownames(x=object)
    dr <- object@reductions

    gzPath <- file.path(dir, "exprMatrix.tsv.gz")
    if ((((ncol(counts)/1000)*(nrow(counts)/1000))>2000) && is(counts, 'sparseMatrix')){
        writeSparseMatrix(counts, gzPath);
    } else {
        mat <- as.matrix(counts)
        df <- as.data.frame(mat, check.names=FALSE)
        df <- data.frame(gene=genes, df, check.names=FALSE)
        z <- gzfile(gzPath, "w")
        write.table(x=df, sep="\t", file=z, quote=FALSE, row.names=FALSE)
        close(con=z)
    }
    embeddings = names(dr)
    embeddings.conf <- c()
    for (embedding in embeddings) {
        emb <- dr[[embedding]]
        df <- emb@cell.embeddings
        if (ncol(df) > 2){
            df <- df[, 1:2]
        }
        colnames(df) <- c("x", "y")
        df <- data.frame(cellId=rownames(df), df, check.names=FALSE)
        write.table(
            df[cellOrder,],
            sep="\t",
            file=file.path(dir, sprintf("%s.coords.tsv", embedding)),
            quote=FALSE,
            row.names=FALSE
        )
        embeddings.conf <- c(
            embeddings.conf,
            sprintf('{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}', embedding)
        )
    }
    df <- data.frame(row.names=cellOrder, check.names=FALSE)
    enum.fields <- c()
    for (i in 1:length(meta.fields)){
        field <- meta.fields[i]
        name <- meta.fields.names[i]
        if (field %in% colnames(meta)){
            df[[name]] <- meta[[field]]
            if (!is.numeric(df[[name]])) {
                enum.fields <- c(enum.fields, name)
            }
        }
    }
    df <- data.frame(Cell=rownames(df), df, check.names=FALSE)
    write.table(
        as.matrix(df[cellOrder,]),
        sep="\t",
        file=file.path(dir, "meta.tsv"),
        quote=FALSE,
        row.names=FALSE
    )
    if (!is.null(features)){
        write.table(
            features,
            sep="\t",
            file=file.path(dir, "quickGenes.tsv"),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE
        )
    }
    config <- '
name="cellbrowser"
shortLabel="cellbrowser"
geneLabel="Feature"
exprMatrix="exprMatrix.tsv.gz"
meta="meta.tsv"
radius=3
alpha=0.5
geneIdType="auto"
clusterField="%s"
labelField="%s"
enumFields=%s
coords=%s'
    enum.string <- paste0("[", paste(paste0('"', enum.fields, '"'), collapse=", "), "]")
    coords.string <- paste0("[", paste(embeddings.conf, collapse = ",\n"), "]")
    config <- sprintf(
        config,
        cluster.field,
        cluster.field,
        enum.string,
        coords.string
    )
    if (!is.null(features)){
        config <- paste(config, 'quickGenesFile="quickGenes.tsv"', sep="\n")
    }
    confPath = file.path(dir, "cellbrowser.conf")
    cat(config, file=confPath)
    cb.dir <- paste(dir, "html_data", sep="/")
    cb <- reticulate::import(module = "cellbrowser")
    cb$cellbrowser$build(dir, cb.dir)
}


export_cellbrowser_data <- function(seurat_data, assay, matrix_slot, resolution, features, rootname){
    tryCatch(
        expr = {
            backup_assay <- DefaultAssay(seurat_data)
            DefaultAssay(seurat_data) <- assay
            meta_fields       <- c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "frip", "blacklisted_fraction")
            meta_fields_names <- c("GEX UMIs",   "Genes",        "Mitochondrial %", "Novelty score",            "ATAC UMIs",   "Peaks",         "Nucleosome signal", "FRiP", "Blacklisted fraction")
            meta_fields <- c(
                meta_fields,
                paste("wsnn_res", resolution, sep=".")
            )
            meta_fields_names <- c(
                meta_fields_names,
                paste("Clustering (", resolution, ")", sep="")
            )
            ExportToCellbrowser(
                seurat_data,
                matrix.slot=matrix_slot,
                dir=rootname,
                cluster.field=paste("Clustering (", resolution, ")", sep="")[1],    # take only the first of possible multiple resolutions
                features=features,
                meta.fields=meta_fields,
                meta.fields.names=meta_fields_names
            )
            print(paste("Export UCSC Cellbrowser data to", rootname, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export UCSC Cellbrowser data to", rootname, sep=" "))
        },
        finally = {
            DefaultAssay(seurat_data) <- backup_assay
        }
    )
}


add_qc_metrics <- function(seurat_data, blacklisted_data, args) {
    backup_assay <- DefaultAssay(seurat_data)

    DefaultAssay(seurat_data) <- "RNA"
    seurat_data$log10_gene_per_log10_umi <- log10(seurat_data$nFeature_RNA) / log10(seurat_data$nCount_RNA)
    seurat_data$mito_percentage <- PercentageFeatureSet(seurat_data, pattern=args$mitopattern) 

    DefaultAssay(seurat_data) <- "ATAC"

    # Nucleosome banding pattern
    #   The histogram of DNA fragment sizes (determined from the paired-end sequencing reads)
    #   should exhibit a strong nucleosome banding pattern corresponding to the length of DNA
    #   wrapped around a single nucleosome. We calculate this per single cell, and quantify
    #   the approximate ratio of mononucleosomal to nucleosome-free fragments (stored as 
    #   nucleosome_signal)
    seurat_data <- NucleosomeSignal(seurat_data, verbose=FALSE)

    # Fraction of fragments in peaks.
    #   Represents the fraction of all fragments that fall within ATAC-seq peaks. Cells with low values
    #   (i.e. <15-20%) often represent low-quality cells or technical artifacts that should be removed.
    #   Note that this value can be sensitive to the set of peaks used.
    fragments_data <- CountFragments(
        fragments=args$fragments,
        cells=colnames(seurat_data),           # limit it to only those cells that are present in seurat_data
        verbose=FALSE
    ) %>% column_to_rownames("CB")             # for easy access to cell barcodes
    seurat_data$fragments <- fragments_data[colnames(seurat_data), "frequency_count"]  # select by rownames to make sure the cells order wasn't accidentally changed
    seurat_data <- FRiP(
        seurat_data,
        assay="ATAC",                          # FRiP can't take the default assay, so we set it explicitly
        total.fragments="fragments",
        col.name="frip",
        verbose=FALSE
    )

    # Ratio reads in genomic blacklist regions.
    #   The ENCODE project has provided a list of blacklist regions, representing reads which are often
    #   associated with artefactual signal. Cells with a high proportion of reads mapping to these areas
    #   often represent technical artifacts and should be removed.
    if (!is.null(blacklisted_data)){
        seurat_data$blacklisted_fraction <- FractionCountsInRegion(
            seurat_data, 
            assay="ATAC",
            regions=blacklisted_data
        )
    } else {
        seurat_data$blacklisted_fraction <- 0  # blacklisted regions file wasn't provided, so we set everything to 0
    }

    print(head(seurat_data@meta.data))
    DefaultAssay(seurat_data) <- backup_assay
    return (seurat_data)
}


get_scaled_norm_seurat_data <- function(seurat_data, args) {
    scaled_norm_seurat_data <- SCTransform(
        seurat_data,
        assay="RNA",
        new.assay.name="SCT",
        variable.features.n=args$highvarcount,
        verbose=FALSE
    )
    return (scaled_norm_seurat_data)
}


apply_qc_filters <- function(seurat_data, args) {
    filtered_seurat_data <- subset(
        seurat_data,
        subset=(nFeature_RNA >= args$mingenes) &
               (nFeature_RNA <= args$maxgenes) &
               (nCount_RNA >= args$gexminumi) &
               (log10_gene_per_log10_umi >= args$minnovelty) &
               (mito_percentage <= args$maxmt) &
               (nCount_ATAC >= args$atacminumi) &
               (nucleosome_signal <= args$maxnuclsignal) &
               (frip >= args$minfrip) &
               (blacklisted_fraction <= args$maxblacklisted)
    )
    return (filtered_seurat_data)
}


export_geom_bar_plot <- function(data, rootname, x_axis, color_by, x_label, y_label, legend_title, plot_title, palette="Paired", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                geom_bar(colour="black") +
                geom_text(stat="count", aes(label=..count..), vjust=-1) +
                xlab(x_label) +
                ylab(y_label) +
                guides(fill=guide_legend(legend_title), x=guide_axis(angle=45)) +
                ggtitle(plot_title) +
                scale_fill_brewer(palette=palette)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export geom bar plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export geom bar plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_geom_density_plot <- function(data, rootname, x_axis, color_by, facet_by, x_left_intercept, x_label, y_label, legend_title, plot_title, x_right_intercept=NULL,  scale_x_log10=FALSE, scale_y_log10=FALSE, zoom_on_intercept=FALSE, alpha=0.9, show_ranked=FALSE, ranked_x_label="Ranked cells", palette="Paired", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            palette_colors <- RColorBrewer::brewer.pal(99, palette)  # use 99 to get all available colors. As we added LETTERS, the order will be correct
            intercept_data <- data %>%
                              dplyr::select(all_of(color_by), all_of(facet_by)) %>%
                              distinct() %>%
                              arrange(all_of(color_by)) %>%
                              add_column(color=palette_colors[1:nrow(.)],x_left=x_left_intercept)

            plot <- ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                    geom_density(alpha=alpha) +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(fill=guide_legend(legend_title)) +
                    ggtitle(plot_title) +
                    facet_wrap(as.formula(paste("~", facet_by))) +
                    scale_fill_brewer(palette=palette) +
                    geom_vline(intercept_data, mapping=aes(xintercept=x_left), color=intercept_data$color, alpha=0.7) +
                    geom_label_repel(
                        intercept_data, mapping=aes(x=x_left, y=Inf, label=x_left),
                        color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                        show.legend=FALSE
                    )

            if (!is.null(x_right_intercept)){
                intercept_data <- intercept_data %>% add_column(x_right=x_right_intercept)
                plot <- plot +
                        geom_vline(intercept_data, mapping=aes(xintercept=x_right), color=intercept_data$color, alpha=0.7) +
                        geom_label_repel(
                            intercept_data, mapping=aes(x=x_right, y=Inf, label=x_right),
                            color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                            show.legend=FALSE
                        )
            }

            if (scale_x_log10){ plot <- plot + scale_x_log10() }
            if (scale_y_log10){ plot <- plot + scale_y_log10() }

            if (zoom_on_intercept) {
                zoomed_plot <- ggplot(data, aes_string(x=x_axis, color=color_by)) +
                               geom_density(show.legend=FALSE, size=2) +
                               xlab(x_label) +
                               ylab(y_label) +
                               scale_color_brewer(palette=palette) +
                               geom_vline(intercept_data, mapping=aes(xintercept=x_left), color=intercept_data$color, alpha=0.7) +
                               geom_label_repel(
                                   intercept_data, mapping=aes(x=x_left, y=Inf, label=x_left),
                                   color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                                   show.legend=FALSE
                               )
                if (!is.null(x_right_intercept)){
                    zoomed_plot <- zoomed_plot +
                                   coord_cartesian(xlim=c(min(intercept_data$x_left), max(intercept_data$x_right))) +
                                   geom_vline(intercept_data, mapping=aes(xintercept=x_right), color=intercept_data$color, alpha=0.7) +
                                   geom_label_repel(
                                       intercept_data, mapping=aes(x=x_right, y=Inf, label=x_right),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                                       show.legend=FALSE
                                   )
                } else {
                    zoomed_plot <- zoomed_plot + coord_cartesian(xlim=c(NA, max(intercept_data$x_left)))
                }
                if (scale_x_log10){ zoomed_plot <- zoomed_plot + scale_x_log10() }
                if (scale_y_log10){ zoomed_plot <- zoomed_plot + scale_y_log10() }
                plot <- plot / zoomed_plot
            }

            if (show_ranked) {
                ranked_plot <- data %>%
                               arrange(get(color_by), get(x_axis)) %>%
                               ggplot(aes_string(x=seq_along(data[[x_axis]]), y=x_axis, color=color_by)) +
                               geom_point(show.legend=FALSE, size=0.5) +
                               xlab(ranked_x_label) +
                               ylab(x_label) +
                               scale_y_log10() +
                               scale_color_brewer(palette=palette) +
                               geom_hline(intercept_data, mapping=aes(yintercept=x_left), color=intercept_data$color, alpha=0.7) +
                               geom_label_repel(
                                   intercept_data, mapping=aes(x=Inf, y=x_left, label=x_left),
                                   color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="x", size=3,
                                   show.legend=FALSE
                               )
                if (!is.null(x_right_intercept)){
                    ranked_plot <- ranked_plot +
                                   geom_hline(intercept_data, mapping=aes(yintercept=x_right), color=intercept_data$color, alpha=0.7) +
                                   geom_label_repel(
                                       intercept_data, mapping=aes(x=-Inf, y=x_right, label=x_right),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="x", size=3,
                                       show.legend=FALSE
                                   )
                }
                plot <- plot / ranked_plot
            }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export geom density plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export geom density plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_geom_point_plot <- function(data, rootname, x_axis, y_axis, facet_by, x_left_intercept, y_low_intercept, color_by, colors, color_limits, color_break, x_label, y_label, legend_title, plot_title, y_high_intercept=NULL, scale_x_log10=FALSE, scale_y_log10=FALSE, alpha=0.2, alpha_intercept=0.5, palette="Paired", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            palette_colors <- RColorBrewer::brewer.pal(99, palette)  # use 99 to get all available colors. As we added LETTERS, the order will be correct
            intercept_data <- data %>%
                              dplyr::select(all_of(facet_by)) %>%
                              distinct() %>%
                              arrange(all_of(facet_by)) %>%
                              add_column(color=palette_colors[1:nrow(.)], x_left=x_left_intercept, y_low=y_low_intercept)

            if (!is.null(y_high_intercept)){
                intercept_data <- intercept_data %>% add_column(y_high=y_high_intercept)
            }

            plot <- ggplot(data, aes_string(x=x_axis, y=y_axis, color=color_by)) +
                    geom_point(alpha=alpha) +
                    scale_colour_gradientn(
                        colours=c(colors[1], colors),
                        values=rescale(c(color_limits[1], color_break-0.01*color_break, color_break, color_limits[2])),
                        breaks=c(color_break),
                        limits=color_limits
                    ) +
                    stat_smooth(method=lm) +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(color=guide_colourbar(legend_title)) +
                    ggtitle(plot_title) +
                    facet_wrap(as.formula(paste("~", facet_by))) +
                    geom_vline(intercept_data, mapping=aes(xintercept=x_left), color=intercept_data$color, alpha=alpha_intercept) +
                    geom_hline(intercept_data, mapping=aes(yintercept=y_low), color=intercept_data$color, alpha=alpha_intercept) +
                    geom_label_repel(
                        intercept_data, mapping=aes(x=x_left, y=Inf, label=x_left),
                        color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="y", size=3,
                        show.legend=FALSE
                    ) +
                    geom_label_repel(
                        intercept_data, mapping=aes(x=Inf, y=y_low, label=y_low),
                        color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="x", size=3,
                        show.legend=FALSE
                    )

            if (!is.null(y_high_intercept)){
                plot <- plot +
                        geom_hline(intercept_data, mapping=aes(yintercept=y_high), color=intercept_data$color, alpha=alpha_intercept) +
                        geom_label_repel(
                            intercept_data, mapping=aes(x=Inf, y=y_high, label=y_high),
                            color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="x", size=3,
                            show.legend=FALSE
                        )
            }

            if (scale_x_log10){ plot <- plot + scale_x_log10() }
            if (scale_y_log10){ plot <- plot + scale_y_log10() }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export geom point plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export geom point plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_vln_plot <- function(data, features, labels, rootname, plot_title, legend_title, log=FALSE, group_by=NULL, hide_x_text=FALSE, pt_size=NULL, palette=NULL, combine_guides=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plots <- VlnPlot(
                         data,
                         features=features,
                         pt.size=pt_size,
                         group.by=group_by,
                         log=log,
                         combine=FALSE       # to return a list of gglots
                     )
            plots <- lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] +
                              ggtitle(labels[i]) +
                              theme_gray() +
                              theme(axis.title.x=element_blank()) +
                              guides(fill=guide_legend(legend_title)) +
                              stat_boxplot(width=0.15, geom="errorbar") +
                              geom_boxplot(width=0.15, outlier.alpha=0) +
                              RotatedAxis()
                if (!is.null(palette)){ plots[[i]] <- plots[[i]] + scale_fill_brewer(palette=palette) }
                if (hide_x_text){ plots[[i]] <- plots[[i]] + theme(axis.text.x=element_blank()) }
                return (plots[[i]])
            })
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export violin plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export violin plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_dim_plot <- function(data, rootname, reduction, plot_title, legend_title, split_by=NULL, group_by=NULL, perc_split_by=NULL, perc_group_by=NULL, label=FALSE, palette=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- DimPlot(
                        data,
                        reduction=reduction,
                        split.by=split_by,
                        group.by=group_by,
                        label=label
                    ) +
                    theme_gray() +
                    ggtitle(plot_title) +
                    guides(color=guide_legend(legend_title, override.aes=list(size=3)))

            if (!is.null(palette)){ plot <- plot + scale_color_brewer(palette=palette) }

            if(!is.null(perc_split_by) && !is.null(perc_group_by)){
                width <- 2 * width
                perc_data <- data@meta.data %>%
                             dplyr::group_by(across(all_of(perc_split_by)), across(all_of(perc_group_by))) %>%
                             summarise(counts=n(), .groups="drop_last") %>%               # drop the perc_group_by grouping level so we can get only groups defined by perc_split_by
                             mutate(freq=counts/sum(counts)*100) %>%                      # sum is taken for the group defined by perc_split_by
                             arrange(all_of(perc_split_by), all_of(perc_group_by))        # sort for consistency
                label_data <- data@meta.data %>%
                              dplyr::group_by(across(all_of(perc_split_by))) %>%
                              summarise(counts=n(), .groups="drop_last") %>%              # drops all grouping as we have only one level
                              arrange(all_of(perc_split_by))                              # sort for consistency
                perc_plot <- ggplot(perc_data, aes_string(x=perc_split_by, y="freq", fill=perc_group_by)) +
                             geom_col(position="dodge", width=0.9, linetype="solid", color="black", show.legend=FALSE) +
                             xlab("") +
                             ylab("Cells percentage") +
                             theme_gray() +
                             geom_label_repel(
                                label_data, mapping=aes(y=-Inf, label=counts),
                                color="black", fill="white", segment.colour=NA,
                                direction="y", size=3, show.legend=FALSE
                             ) +
                             RotatedAxis()
                plot <- plot + perc_plot
            }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export dim plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export dim plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_elbow_plot <- function(data, rootname, plot_title, ndims=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- ElbowPlot(
                        data,
                        ndims=ndims,
                        reduction="pca"
                    ) +
                    theme_gray() +
                    ggtitle(plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export elbow plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export elbow plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_tss_plot <- function(data, rootname, plot_title, group_by_value=NULL, idents=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            group_by <- NULL
            if (!is.null(group_by_value)){
                data$tss_group_by <- ifelse(
                    data$TSS.enrichment >= group_by_value ,
                    paste("High ", "(bigger or equal to ", group_by_value, ")", sep=""),
                    paste("Low ", "(smaller than ", group_by_value, ")", sep="")
                )
                group_by <- "tss_group_by"
            }
            plot <- TSSPlot(
                        data,
                        group.by=group_by,
                        idents=idents
                    ) +
                    theme_gray() +
                    ggtitle(plot_title) +
                    NoLegend()

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export TSS Enrichment plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export TSS Enrichment plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_fragments_hist <- function(data, rootname, plot_title, group_by_value=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            group_by <- NULL
            if (!is.null(group_by_value)){
                data$ns_group_by <- ifelse(
                    data$nucleosome_signal >= group_by_value ,
                    paste("High nucleosome signal ", "(bigger or equal to ", group_by_value, ")", sep=""),
                    paste("Low nucleosome signal", "(smaller than ", group_by_value, ")", sep="")
                )
                group_by <- "ns_group_by"
            }
            plot <- FragmentHistogram(
                        data,
                        group.by=group_by
                    ) +
                    theme_gray() +
                    ggtitle(plot_title) +
                    NoLegend()

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export fragments length histogram to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export fragments length histogram to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_all_qc_plots <- function(seurat_data, suffix, args){
    export_geom_bar_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "cell_count", sep="_"),
        x_axis="orig.ident",
        color_by="orig.ident",
        x_label="Identity",
        y_label="Cells",
        legend_title="Identity",
        plot_title=paste("Number of cells per dataset (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gex_umi_dnst", sep="_"),
        x_axis="nCount_RNA",
        color_by="orig.ident",
        facet_by="orig.ident",
        x_left_intercept=args$gexminumi,
        x_label="UMIs per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("GEX UMI density per cell (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "atac_umi_dnst", sep="_"),
        x_axis="nCount_ATAC",
        color_by="orig.ident",
        facet_by="orig.ident",
        x_left_intercept=args$atacminumi,
        x_label="UMIs per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("ATAC UMI density per cell (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gene_dnst", sep="_"),
        x_axis="nFeature_RNA",
        color_by="orig.ident",
        facet_by="orig.ident",
        x_left_intercept=args$mingenes,
        x_right_intercept=args$maxgenes,
        x_label="Genes per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Gene density per cell (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "peak_dnst", sep="_"),
        x_axis="nFeature_ATAC",
        color_by="orig.ident",
        facet_by="orig.ident",
        x_left_intercept=0,
        # x_right_intercept=args$maxgenes,
        x_label="Peaks per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Peak density per cell (", suffix, ")", sep=""),
        scale_x_log10=FALSE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "bl_cnts_dnst", sep="_"),
        x_axis="blacklisted_fraction",
        color_by="orig.ident",
        facet_by="orig.ident",
        x_left_intercept=args$maxblacklisted,
        x_label="Fraction of reads within blacklisted regions per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Density of fraction of reads within blacklisted regions per cell (", suffix, ")", sep=""),
        scale_x_log10=FALSE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        pdf=args$pdf
    )
    export_geom_point_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gex_atac_umi_corr", sep="_"),
        facet_by="orig.ident",
        x_axis="nCount_ATAC",
        x_label="ATAC UMIs per cell",
        y_axis="nCount_RNA",
        y_label="GEX UMIs per cell",
        x_left_intercept=args$atacminumi,
        y_low_intercept=args$gexminumi,
        alpha_intercept=1,
        color_by="mito_percentage",
        colors=c("lightslateblue", "red", "green"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        legend_title="Mitochondrial %",
        plot_title=paste("GEX vs ATAC UMIs per cell correlation (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        pdf=args$pdf
    )
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "ATAC"
    export_fragments_hist(
        data=seurat_data,
        group_by_value=args$maxnuclsignal,
        rootname=paste(args$output, suffix, "frg_len_hist", sep="_"),
        plot_title=paste("Fragments Length Histogram (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    DefaultAssay(seurat_data) <- backup_assay
    export_geom_point_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gene_umi_corr", sep="_"),
        x_axis="nCount_RNA",
        y_axis="nFeature_RNA",
        facet_by="orig.ident",
        x_left_intercept=args$gexminumi,
        y_low_intercept=args$mingenes,
        y_high_intercept=args$maxgenes,
        color_by="mito_percentage",
        colors=c("lightslateblue", "red", "green"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        x_label="UMIs per cell",
        y_label="Genes per cell",
        legend_title="Mitochondrial %",
        plot_title=paste("Genes vs GEX UMIs per cell correlation (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "mito_perc_dnst", sep="_"),
        x_axis="mito_percentage",
        color_by="orig.ident",
        facet_by="orig.ident",
        x_left_intercept=args$maxmt,
        x_label="Percentage of transcripts mapped to mitochondrial genes per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Density of transcripts mapped to mitochondrial genes per cell (", suffix, ")", sep=""),
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "nvlt_score_dnst", sep="_"),
        x_axis="log10_gene_per_log10_umi",
        color_by="orig.ident",
        facet_by="orig.ident",
        x_left_intercept=args$minnovelty,
        x_label="log10 Genes / log10 UMIs per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Novelty score density per cell (", suffix, ")", sep=""),
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    export_vln_plot(
        data=seurat_data,
        features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "frip", "blacklisted_fraction"),
        labels=c(  "GEX UMIs",   "Genes",        "Mitochondrial %", "Novelty score",            "ATAC UMIs",   "Peaks",         "Nucleosome signal", "FRiP", "Fr. of reads in bl-ted reg."),
        rootname=paste(args$output, suffix, "qc_mtrcs", sep="_"),
        plot_title=paste("QC metrics densities per cell (", suffix, ")", sep=""),
        legend_title="Identity",
        hide_x_text=TRUE,
        pt_size=0,
        palette="Paired",
        combine_guides="collect",
        pdf=args$pdf
    )
}


export_all_clustering_plots <- function(seurat_data, suffix, args) {
    cluster_prefix <- "wsnn"
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        export_dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title=paste("Clustered UMAP projected PCA of filtered GEX datasets. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "gex_umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="atacumap",
            plot_title=paste("Clustered UMAP projected LSI of filtered ATAC datasets. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "atac_umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="wnnumap",
            plot_title=paste("Clustered UMAP projected WNN. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "wnn_umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
    }
}


export_all_dimensionality_plots <- function(seurat_data, suffix, args) {
    export_elbow_plot(
        data=seurat_data,
        ndims=50,
        rootname=paste(args$output, suffix, "elbow", sep="_"),
        plot_title="Elbow plot from PCA of filtered integrated/scaled datasets",
        pdf=args$pdf
    )
    export_dim_plot(
        data=seurat_data,
        reduction="pca",
        plot_title="PCA of filtered integrated/scaled datasets",
        legend_title="Identity",
        rootname=paste(args$output, suffix, "pca", sep="_"),
        palette="Paired",
        pdf=args$pdf
    )
}


export_dot_plot <- function(data, features, rootname, plot_title, x_label, y_label, cluster_idents=FALSE, min_pct=0.01, col_min=-2.5, col_max=2.5, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- DotPlot(
                        data,
                        features=features,
                        cluster.idents=cluster_idents,
                        dot.min=min_pct,
                        col.min=col_min,
                        col.max=col_max,
                        scale=TRUE,
                        scale.by="size"  # for optimal perception
                    ) +
                    xlab(x_label) +
                    ylab(y_label) +
                    theme_gray() +
                    ggtitle(plot_title) +
                    RotatedAxis()

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export dot plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export dot plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_feature_plot <- function(data, features, labels, rootname, reduction, plot_title, split_by=NULL, label=FALSE, order=FALSE, min_cutoff=NA, max_cutoff=NA, pt_size=NULL, combine_guides=NULL, alpha=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plots <- FeaturePlot(
                        data,
                        features=features,
                        pt.size=pt_size,
                        order=order,
                        min.cutoff=min_cutoff,
                        max.cutoff=max_cutoff,
                        reduction=reduction,
                        split.by=split_by,
                        label=label,
                        combine=FALSE       # to return a list of gglots
                    )
            plots <- lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] + ggtitle(labels[i]) + theme_gray()
                if (!is.null(alpha)) { plots[[i]]$layers[[1]]$aes_params$alpha <- alpha }
                return (plots[[i]])
            })
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export feature plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export feature plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_all_expression_plots <- function(seurat_data, suffix, args, assay="RNA") {
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        Idents(seurat_data) <- paste("wsnn_res", current_resolution, sep=".")
        export_dot_plot(
            data=seurat_data,
            features=args$gexfeatures,
            plot_title=paste("Scaled average log normalized gene expression per cluster. Resolution", current_resolution),
            x_label="Genes",
            y_label="Clusters",
            cluster_idents=FALSE,
            rootname=paste(args$output, suffix, "avg_per_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_feature_plot(
            data=seurat_data,
            features=args$gexfeatures,
            labels=args$gexfeatures,
            reduction="wnnumap",
            plot_title=paste("Log normalized gene expression per cell of clustered datasets. Resolution", current_resolution),
            label=TRUE,
            order=TRUE,
            max_cutoff="q99",  # to prevent cells with overexpressed gene from distoring the color bar
            rootname=paste(args$output, suffix, "per_clst_cell_res", current_resolution, sep="_"),
            combine_guides="keep",
            pdf=args$pdf
        )
        export_vln_plot(
            data=seurat_data,
            features=args$gexfeatures,
            labels=args$gexfeatures,
            plot_title=paste("Log normalized gene expression densities per cluster. Resolution", current_resolution),
            legend_title="Cluster",
            log=TRUE,
            pt_size=0,
            combine_guides="collect",
            rootname=paste(args$output, suffix, "dnst_per_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        Idents(seurat_data) <- "orig.ident"
    }
    DefaultAssay(seurat_data) <- backup_assay
}


export_rds <- function(data, location){
    tryCatch(
        expr = {
            saveRDS(data, location)
            print(paste("Export data as RDS to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data as RDS to", location, sep=" "))
        }
    )
}


load_seurat_data <- function(args) {
    raw_data <- Read10X(data.dir=args$mex) 
    seurat_data <- CreateSeuratObject(
        counts=raw_data$`Gene Expression`,
        min.cells=args$gexmincells,
        names.delim="-",
        names.field=2
    )

    annotation <- rtracklayer::import(args$annotations, format="GFF")

    seurat_data[["ATAC"]] <- CreateChromatinAssay(
        counts=raw_data$Peaks,
        sep=c(":", "-"),
        fragments=args$fragments,
        min.cells=args$atacmincells,
        annotation=annotation
    )

    if (args$callpeaks){
        print("Replacing Cell Ranger's peaks with MACS2 peaks.")
        backup_assay <- DefaultAssay(seurat_data)
        DefaultAssay(seurat_data) <- "ATAC"
        # We use group.by="orig.ident" to force CallPeaks using only those fragments
        # that belong to the cells already identified by Cell Ranger as valid.
        # Otherwise, CallPeaks function will use all fragments even from invalid cells.
        # Grouping by orig.ident will results in only 1 group, which is fine as long
        # as our input mex matrix includes only one dataset (not aggregated).
        macs2_peaks <- CallPeaks(seurat_data, group.by="orig.ident")
        macs2_counts <- FeatureMatrix(
            fragments=Fragments(seurat_data),
            sep=c(":", "-"),
            features=macs2_peaks,
            cells=colnames(seurat_data),
            verbose=FALSE,
        )
        atac_assay <- CreateChromatinAssay(
            counts=macs2_counts,
            sep=c(":", "-"),
            fragments=Fragments(seurat_data),
            min.cells=args$atacmincells,
            min.features=-1,              # as they check ncount.cell > min.features and by default it's 0, we will remove cells without peaks and won't be able to add new assay to our seurat_data
            annotation=annotation
        )
        seurat_data[["ATAC"]] <- atac_assay
        DefaultAssay(seurat_data) <- backup_assay
    }

    Idents(seurat_data) <- "orig.ident"
    return (seurat_data)
}


load_blacklisted_data <- function(location) {
    default_blacklisted_data <- NULL
    if (!is.null(location)){
        blacklisted_data <- rtracklayer::import(location, format="BED")
        print(paste("Blacklisted regions data is successfully loaded from ", location))
        return (blacklisted_data)
    }
    print("Blacklisted regions data is not provided. Cells won't be filtered by --maxblacklisted")
    return (default_blacklisted_data)
}


get_args <- function(){
    parser <- ArgumentParser(description='Runs Seurat for comparative scRNA-seq analysis of across experimental conditions')
    # Loading data
    parser$add_argument(
        "--mex",
        help=paste(
            "Path to the folder with feature-barcode matrix from Cell Ranger ARC Count",
            "in MEX format. The rows consist of all the gene and peak features concatenated",
            "together and the columns are restricted to those barcodes that are identified",
            "as cells."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment observed in",
            "the experiment in TSV format. Tbi-index file is required."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--annotations",
        help="Path to the genome annotation file in GTF format",
        type="character", required="True"
    )
    parser$add_argument(
        "--blacklisted",
        help="Path to the blacklisted regions file in BED format",
        type="character"
    )
    # Filtering
    parser$add_argument(
        "--gexmincells",
        help=paste(
            "Include only GEX features detected in at least this many cells.",
            "Default: 5"
        ),
        type="integer", default=5
    )
    parser$add_argument(
        "--mingenes",
        help=paste(
            "Include cells where at least this many GEX features are detected",
            "Default: 250"
        ),
        type="integer", default=250
    )
    parser$add_argument(
        "--maxgenes",
        help=paste(
            "Include cells with the number of GEX features not bigger than this value.",
            "Default: 5000"
        ),
        type="integer", default=5000
    )
    parser$add_argument(
        "--gexminumi",
        help=paste(
            "Include cells where at least this many GEX UMIs (transcripts) are detected.",
            "Default: 500"
        ),
        type="integer", default=500
    )
    parser$add_argument(
        "--mitopattern",
        help="Regex pattern to identify mitochondrial GEX features. Default: '^Mt-'",
        type="character", default="^Mt-"
    )
    parser$add_argument(
        "--maxmt",
        help=paste(
            "Include cells with the percentage of GEX transcripts mapped to mitochondrial",
            "genes not bigger than this value. Default: 5"
        ),
        type="double", default=5
    )
    parser$add_argument(
        "--minnovelty",
        help=paste(
            "Include cells with the novelty score not lower than this value,",
            "calculated for GEX as log10(genes)/log10(UMIs).",
            "Default: 0.8"
        ),
        type="double", default=0.8
    )

    parser$add_argument(
        "--atacmincells",
        help=paste(
            "Include only ATAC features detected in at least this many cells.",
            "Default: 5"
        ),
        type="integer", default=5
    )
    parser$add_argument(
        "--atacminumi",
        help=paste(
            "Include cells where at least this many ATAC UMIs (transcripts) are detected.",
            "Default: 1000"
        ),
        type="integer", default=1000
    )
    parser$add_argument(
        "--maxnuclsignal",
        help=paste(
            "Include cells with the nucleosome signal not bigger than this value.",
            "Nucleosome signal quantifies the approximate ratio of mononucleosomal",
            "to nucleosome-free fragments.",
            "Default: 4"
        ),
        type="double", default=4
    )
    parser$add_argument(
        "--minfrip",
        help=paste(
            "Include cells with the FRiP not lower than this value.",
            "Default: 0.15"
        ),
        type="double", default=0.15
    )
    parser$add_argument(
        "--maxblacklisted",
        help=paste(
            "Include cells with the ratio of reads in genomic blacklist regions",
            "not bigger than this value.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--callpeaks",
        help=paste(
            "Call peaks with MACS2 instead of those that are provided by Cell Ranger ARC Count.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--gexfeatures",
        help="GEX features of interest to evaluate expression. Default: None",
        type="character", nargs="*"
    )
    parser$add_argument(
        "--highvarcount",
        help=paste(
            "Number of highly variable features to detect. Used for datasets integration,",
            "scaling, and dimensional reduction. Default: 3000"
        ),
        type="integer", default=3000
    )
    parser$add_argument(
        "--gexndim",
        help=paste(
            "Number of principal components to use in GEX UMAP projection and clustering",
            "(from 1 to 50). Use Elbow plot to adjust this parameter. Default: 50"
        ),
        type="integer", default=50
    )
    parser$add_argument(
        "--atacndim",
        help=paste(
            "Number of principal components to use in ATAC UMAP projection and clustering",
            "(from 1 to 50). Use Elbow plot to adjust this parameter. Default: 50"
        ),
        type="integer", default=50
    )
    parser$add_argument(
        "--resolution",
        help="Clustering resolution. Can be set as an array. Default: 0.3",
        type="double", default=c(0.3), nargs="*"
    )
    # Export results
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--rds",
        help="Save Seurat data to RDS file. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./seurat",
        type="character", default="./seurat"
    )

    # Performance parameters
    parser$add_argument(
        "--threads",
        help="Threads. Default: 1",
        type="integer", default=1
    )

    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}

args <- get_args()


cat("Step 0: Runtime configuration\n")


print(args)
print(paste("Setting parallelizations threads to", args$threads))
set_threads(args$threads)


cat("\n\nStep 1: Loading raw datasets\n")
print(paste("Loading gene/peak-barcode matrices from", args$mex))
print(paste("Loading fragments from", args$fragments))
print(paste("Loading annotations from from", args$annotations))
seurat_data <- load_seurat_data(args)
blacklisted_data <- load_blacklisted_data(args$blacklisted)
print("Adding QC metrics to raw seurat data")
seurat_data <- add_qc_metrics(seurat_data, blacklisted_data, args)
export_all_qc_plots(seurat_data, "raw", args)

cat("\n\nStep 2: Filtering raw datasets\n")
seurat_data <- apply_qc_filters(seurat_data, args)
export_all_qc_plots(seurat_data, "fltr", args)

cat("\n\nStep 3: Running GEX analysis\n")
backup_assay <- DefaultAssay(seurat_data)
DefaultAssay(seurat_data) <- "RNA"
seurat_data <- get_scaled_norm_seurat_data(seurat_data, args)
seurat_data <- RunPCA(seurat_data, npcs=50, verbose=TRUE)
seurat_data <- RunUMAP(
    seurat_data,
    reduction="pca",
    dims=1:args$gexndim,
    reduction.name="rnaumap",
    reduction.key="RNAUMAP_",
    verbose=FALSE
)
export_all_dimensionality_plots(seurat_data, "ntgr", args)
DefaultAssay(seurat_data) <- backup_assay

cat("\n\nStep 4: Running ATAC analysis\n")
backup_assay <- DefaultAssay(seurat_data)
DefaultAssay(seurat_data) <- "ATAC"
seurat_data <- RunTFIDF(seurat_data, verbose=FALSE)                           # normalizes across cells to correct for differences in cellular sequencing depth, and across peaks to give higher values to more rare peaks
seurat_data <- FindTopFeatures(seurat_data, min.cutoff="q75", verbose=FALSE)
seurat_data <- RunSVD(seurat_data, verbose=FALSE)
seurat_data <- RunUMAP(
    seurat_data,
    reduction="lsi",
    dims=2:args$atacndim,
    reduction.name="atacumap",
    reduction.key="ATACUMAP_",
    verbose=FALSE
)
DefaultAssay(seurat_data) <- backup_assay

cat("\n\nStep 6: Running WNN analysis\n")
seurat_data <- FindMultiModalNeighbors(
    seurat_data,
    reduction.list=list("pca", "lsi"),
    dims.list=list(1:args$gexndim, 2:args$atacndim),
    snn.graph.name="wsnn",
    weighted.nn.name="weighted.nn",
    verbose=FALSE
)
seurat_data <- RunUMAP(
    seurat_data,
    nn.name="weighted.nn",
    reduction.name="wnnumap",
    reduction.key="WNNUMAP_",
    verbose=FALSE
)
seurat_data <- FindClusters(
    seurat_data,
    graph.name="wsnn",
    algorithm=3,
    resolution=args$resolution,
    verbose=FALSE
)

export_all_clustering_plots(seurat_data, "clst", args)

if (args$rds){ export_rds(seurat_data, paste(args$output, "_clst_data.rds", sep="")) }

if (!is.null(args$gexfeatures)){
    print("Check genes of interest to include only those that are present in the datasets")
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "RNA"
    args$gexfeatures <- unique(args$gexfeatures)
    args$gexfeatures <- args$gexfeatures[args$gexfeatures %in% as.vector(as.character(rownames(seurat_data)))]
    print(args$gexfeatures)
    DefaultAssay(seurat_data) <- backup_assay
}

if (!is.null(args$gexfeatures)){ export_all_expression_plots(seurat_data, "expr", args, assay="RNA") }

print("Exporting UCSC Cellbrowser data")
export_cellbrowser_data(
    seurat_data, assay="RNA", matrix_slot="data",
    resolution=args$resolution, features=args$gexfeatures,
    rootname=paste(args$output, "_cellbrowser", sep="")
)