#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB should be good by default

suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(DESeq2))
suppressMessages(library(tibble))
suppressMessages(library(ggplot2))
suppressMessages(library(garnett))
suppressMessages(library(argparse))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(reticulate))
suppressMessages(library(sctransform))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))


####################################################################################
# For details and treshold values look at
#    https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html
#    https://github.com/satijalab/seurat/issues/1679
#    https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html
#    https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/cellbrowser.html
#    https://github.com/satijalab/seurat-wrappers
####################################################################################

# To use Garnett classifiers trained with Ensembl gene names
SPECIES_DATA <- list(
    "hs"="org.Hs.eg.db",
    "mm"="org.Mm.eg.db"
)

###############################################################################################################
# Modified code from
# https://github.com/maximilianh/cellBrowser/blob/master/src/cbPyLib/cellbrowser/R/ExportToCellbrowser-seurat.R
###############################################################################################################


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
        df[[name]] <- meta[[field]]
        if (!is.numeric(df[[name]])) {
            enum.fields <- c(enum.fields, name)
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
    config <- '
name="cellbrowser"
shortLabel="cellbrowser"
exprMatrix="exprMatrix.tsv.gz"
meta="meta.tsv"
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
    confPath = file.path(dir, "cellbrowser.conf")
    cat(config, file=confPath)
    cb.dir <- paste(dir, "html_data", sep="/")
    cb <- import(module = "cellbrowser")
    cb$cellbrowser$build(dir, cb.dir)
}


###############################################################################################################


set_threads <- function (threads) {
    invisible(capture.output(plan("multiprocess", workers=threads)))
    invisible(capture.output(plan()))
    invisible(capture.output(setDTthreads(threads)))
}


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


load_cell_cycle_data <- function (location) {
    cell_cycle_data <- read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    return (cell_cycle_data)
}


load_cell_identity_data <- function (location) {
    cell_identity_data <- read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    return (cell_identity_data)
}


load_condition_data <- function(location, cell_identity_data) {
    default_condition_data <- data.frame(
        library_id=cell_identity_data$library_id,
        condition=rownames(cell_identity_data),
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    if (!is.null(location)){
        condition_data <- read.table(
            location,
            sep=get_file_type(location),
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )
        if (all(is.element(cell_identity_data$library_id, condition_data$library_id))){
            print(paste("Condition data is successfully loaded from ", location))
            return (condition_data)
        } else {
            print(paste("Applying defaults - failed to load condition data from ", location))
            return (default_condition_data)
        }
    }
    print("Condition data is not provided. Applying defaults")
    return (default_condition_data)
}


load_barcodes_data <- function(location, seurat_data) {
    default_barcodes_data <- Cells(seurat_data)                               # to include all available cells
    if (!is.null(location)){
        barcodes_data <- read.table(
            location,
            sep=get_file_type(location),
            header=FALSE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )[,1]                                                                 # to get it as vector
        print(paste("Barcodes data is successfully loaded from ", location))
        return (barcodes_data)
    }
    print("Barcodes data is not provided. Using all cells")
    return (default_barcodes_data)
}


load_seurat_data <- function(location, mincells, cell_identity_data, condition_data) {
    seurat_data <- CreateSeuratObject(
        counts=Read10X(data.dir=location),
        min.cells=mincells,
        names.delim="-",  # to get cell identity index from Cellranger aggr output
        names.field=2
    )
    seurat_data[["new.ident"]] <- cell_identity_data$library_id[Idents(seurat_data)]
    seurat_data[["condition"]] <- condition_data$condition[match(seurat_data$new.ident, condition_data$library_id)]
    Idents(seurat_data) <- "new.ident"
    return (seurat_data)
}


apply_cell_filters <- function(seurat_data, barcodes_data) {
    filtered_seurat_data <- subset(seurat_data, cells=barcodes_data)
    return (filtered_seurat_data)
}


add_qc_metrics <- function(seurat_data, mitopattern) {
    # Number of genes detected per UMI: this metric with give us an idea of the
    # complexity of our dataset (more genes detected per UMI, more complex our data)
    # Mitochondrial percentage: this metric will give us a percentage of cell reads
    # originating from the mitochondrial genes
    seurat_data$log10_gene_per_log10_umi <- log10(seurat_data$nFeature_RNA) / log10(seurat_data$nCount_RNA)
    seurat_data$mito_percentage <- PercentageFeatureSet(seurat_data, pattern=mitopattern) 
    return (seurat_data)
}


apply_qc_filters <- function(seurat_data, minfeatures, minumi, minnovelty, maxmt) {
    filtered_seurat_data <- subset(
        seurat_data,
        subset = (nFeature_RNA >= minfeatures) & 
                 (nCount_RNA >= minumi) &
                 (log10_gene_per_log10_umi >= minnovelty) &
                 (mito_percentage <= maxmt)
    )
    return (filtered_seurat_data)
}


explore_unwanted_variation <- function(seurat_data, cell_cycle_data, args) {
     temp_seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
     temp_seurat_data <- CellCycleScoring(
         temp_seurat_data,
         s.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="s", "gene_id"]),
         g2m.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="g2/m", "gene_id"])
    )
    mito_quartiles <- quantile(temp_seurat_data@meta.data$mito_percentage, c(0.25, 0.5, 0.75))
    temp_seurat_data@meta.data$mito_factor <- cut(
        temp_seurat_data@meta.data$mito_percentage, 
        breaks=c(-Inf, mito_quartiles[1], mito_quartiles[2], mito_quartiles[3], Inf), 
        labels=c("Low", "Medium", "Medium high", "High")
    )
    temp_seurat_data <- FindVariableFeatures(temp_seurat_data, verbose=FALSE)
    temp_seurat_data <- ScaleData(temp_seurat_data, verbose=FALSE)
    temp_seurat_data <- RunPCA(temp_seurat_data, verbose=FALSE)
    suppressWarnings(
        temp_seurat_data <- RunUMAP(
            temp_seurat_data,
            reduction="pca",
            dims=1:args$ndim,
            verbose=FALSE
        )
    )
    export_dim_plot(
        data=temp_seurat_data,
        reduction="pca",
        split_by="Phase",
        group_by="Phase",
        rootname=paste(args$output, "_filt_unint_cell_cycle_eff_pca_plot", sep="")
    )
    export_dim_plot(
        data=temp_seurat_data,
        reduction="pca",
        split_by="mito_factor",
        group_by="mito_factor",
        rootname=paste(args$output, "_filt_unint_mito_perc_eff_pca_plot", sep="")
    )
    export_dim_plot(
        data=temp_seurat_data,
        reduction="umap",
        split_by="new.ident",
        rootname=paste(args$output, "_filt_unint_umap_plot_spl_by_ident", sep="")
    )
}


integrate_seurat_data <- function(seurat_data, cell_cycle_data, args) {
    splitted_seurat_data <- SplitObject(seurat_data, split.by="new.ident")
    for (i in 1:length(splitted_seurat_data)) {
        vars_to_regress <- NULL
        if (args$regressmt){
            vars_to_regress <- c("mito_percentage")
        }
        splitted_seurat_data[[i]] <- SCTransform(
            splitted_seurat_data[[i]],
            assay="RNA",
            new.assay.name="SCT",
            variable.features.n=args$highvarcount,
            vars.to.regress=vars_to_regress,
            verbose=FALSE
        )
        tryCatch(
            expr = {
                splitted_seurat_data[[i]] <- CellCycleScoring(
                    splitted_seurat_data[[i]],
                    s.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="s", "gene_id"]),
                    g2m.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="g2/m", "gene_id"]),
                    assay="SCT",
                    verbose=FALSE
                )
            },
            error = function(e){
                print("Failed to run cell cycle scoring on SCT normalized data")
            }
        )
        if (args$regresscellcycle){
            if (!is.null(vars_to_regress)) {
                vars_to_regress <- append(vars_to_regress, c("S.Score", "G2M.Score"))
            } else {
                vars_to_regress <- c("S.Score", "G2M.Score")
            }
            splitted_seurat_data[[i]] <- SCTransform(
                splitted_seurat_data[[i]],
                assay="RNA",
                new.assay.name="SCT",
                variable.features.n=args$highvarcount,
                vars.to.regress=vars_to_regress,
                verbose=FALSE
            )            
        }
    }
    set_threads(1)  # need to use one thread, otherwise gets stuck
    integration_features <- SelectIntegrationFeatures(splitted_seurat_data, nfeatures=args$highvarcount)
    splitted_seurat_data <- PrepSCTIntegration(
        splitted_seurat_data, 
        anchor.features=integration_features,
        verbose=FALSE
    )
    integration_anchors <- FindIntegrationAnchors(
        splitted_seurat_data,
        normalization.method="SCT",
        anchor.features=integration_features,
        verbose=FALSE
    )
    integrated_seurat_data <- IntegrateData(
        integration_anchors, 
        new.assay.name="integrated",
        normalization.method="SCT",
        k.weight=min(min(table(Idents(seurat_data))), 100),  # 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering
        verbose=FALSE
    )
    DefaultAssay(integrated_seurat_data) <- "integrated"
    set_threads(args$threads)
    return (integrated_seurat_data)
}


get_all_putative_markers <- function(seurat_data, args, assay="RNA", min_diff_pct=-Inf){
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    all_putative_markers <- NULL
    for (i in 1:length(args$resolution)) {
        resolution <- args$resolution[i]
        Idents(seurat_data) <- paste("integrated_snn_res", resolution, sep=".")
        markers <- FindAllMarkers(
            seurat_data,
            logfc.threshold=args$logfc,
            min.pct=args$minpct,
            only.pos=args$onlypos,
            test.use=args$testuse,
            min.diff.pct=min_diff_pct,
            verbose=FALSE
        ) %>% relocate(cluster, gene, .before=1)
        if (nrow(markers) > 0) {
            markers <- markers %>% cbind(resolution=resolution, .)
        } else {
            markers <- markers %>% add_column(resolution=numeric(), .before=1)  # safety measure in case markers was empty
        }
        if (!is.null(all_putative_markers)) {
            all_putative_markers <- rbind(all_putative_markers, markers)
        } else {
            all_putative_markers <- markers
        }
        Idents(seurat_data) <- "new.ident"
    }
    DefaultAssay(seurat_data) <- backup_assay
    return (all_putative_markers)
}


get_all_conserved_markers <- function(seurat_data, args, assay="RNA", min_diff_pct=-Inf){
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "RNA"
    all_conserved_markers <- NULL
    for (i in 1:length(args$resolution)) {
        resolution <- args$resolution[i]
        clustering_by <- paste("integrated_snn_res", resolution, sep=".")
        Idents(seurat_data) <- clustering_by
        conserved_markers <- map_dfr(
            sort(unique(seurat_data@meta.data[, clustering_by])),
            get_conserved_markers,
            seurat_data,
            "condition",
            resolution,
            args$onlypos,
            args$logfc,
            args$minpct,
            args$testuse,
            min_diff_pct
        )
        if (!is.null(all_conserved_markers)) {
            all_conserved_markers <- rbind(all_conserved_markers, conserved_markers)
        } else {
            all_conserved_markers <- conserved_markers
        }
        Idents(seurat_data) <- "new.ident"
    }
    DefaultAssay(seurat_data) <- backup_assay
    return (all_conserved_markers)
}


get_conserved_markers <- function(cluster, seurat_data, grouping_var, resolution, only_pos=FALSE, logfc_threshold=0.25, min_pct=0.1, test_use="wilcox", min_diff_pct=-Inf){
    conserved_markers <- FindConservedMarkers(
        seurat_data,
        ident.1=cluster,
        grouping.var=grouping_var,
        min.cells.group=0,                 # to catch situation when after splitting by grouping.var we got to few cells in a group
        only.pos=only_pos,
        logfc.threshold=logfc_threshold,
        min.pct=min_pct,
        test.use=test_use,
        min.diff.pct=min_diff_pct,
        verbose=FALSE
    ) %>% rownames_to_column(var="gene")
    if (nrow(conserved_markers) > 0) {
        conserved_markers <- conserved_markers %>% cbind(resolution=resolution, cluster=cluster, .) 
    } else {
        conserved_markers <- conserved_markers %>% add_column(resolution=numeric(), cluster=factor(), .before=1)
    }
    return (conserved_markers)
}


assign_cell_types <- function (seurat_data, classifier, args){
    for (i in 1:length(args$resolution)) {
        resolution <- args$resolution[i]
        monocle_data <- get_monocle_data(
            seurat_data,
            features=VariableFeatures(seurat_data, assay="integrated"),
            cluster_field=paste("integrated_snn_res", resolution, sep=".")
        )
        monocle_data <- classify_cells(
            monocle_data,
            classifier,
            db=get(as.character(SPECIES_DATA[args$species])),
            cluster_extend=TRUE,
            cds_gene_id_type="SYMBOL"
        )
        seurat_data@meta.data[paste("cluster_ext_type_res", resolution, sep=".")] <- monocle_data$cluster_ext_type
    }
    return (seurat_data)
}


get_monocle_data <- function (seurat_data, features, cluster_field, assay="RNA", matrix_slot="data") {
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    gene_metadata <- data.frame(
        "gene_short_name"=features,
        row.names=features,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    cell_metadata <- data.frame(
        seurat_data@meta.data,
        row.names=colnames(seurat_data),
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    cell_metadata["garnett_cluster"] <- cell_metadata[cluster_field]
    expression_matrix <- GetAssayData(seurat_data, slot=matrix_slot)[features, ]
    monocle_data <- new_cell_data_set(
        expression_matrix,
        cell_metadata=cell_metadata,
        gene_metadata=gene_metadata
    )
    DefaultAssay(seurat_data) <- backup_assay
    return (monocle_data)
}


export_geom_bar_plot <- function(data, rootname, x_axis, color_by, x_label, y_label, legend_title, plot_title, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(
                print(
                    ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                    geom_bar() +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(fill=guide_legend(legend_title)) +
                    ggtitle(plot_title)
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(
                    print(
                        ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                        geom_bar() +
                        xlab(x_label) +
                        ylab(y_label) +
                        guides(fill=guide_legend(legend_title)) +
                        ggtitle(plot_title)
                    )
                )
                dev.off()
            }
            print(paste("Export bar plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export bar plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_log10_geom_density_plot <- function(data, rootname, x_axis, color_by, x_intercept, x_label, y_label, legend_title, plot_title, alpha=0.2, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(
                print(
                    ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                    geom_density(alpha=alpha) +
                    scale_x_log10() +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(fill=guide_legend(legend_title)) +
                    geom_vline(xintercept=x_intercept, color="red") +
                    ggtitle(plot_title)
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(
                    print(
                        ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                        geom_density(alpha=alpha) +
                        scale_x_log10() +
                        xlab(x_label) +
                        ylab(y_label) +
                        guides(fill=guide_legend(legend_title)) +
                        geom_vline(xintercept=x_intercept, color="red") +
                        ggtitle(plot_title)
                    )
                )
                dev.off()
            }
            print(paste("Export log10 density plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export log10 density plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_log10_geom_point_plot <- function(data, rootname, x_axis, y_axis, color_by, facet_by, x_intercept, y_intercept, x_label, y_label, legend_title, plot_title, alpha=0.2, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(
                print(
                    ggplot(data, aes_string(x=x_axis, y=y_axis, color=color_by)) +
                    geom_point(alpha=alpha) +
                    scale_colour_gradient(low="gray90", high="black") +
                    stat_smooth(method=lm) +
                    scale_x_log10() + 
                    scale_y_log10() + 
                    geom_vline(xintercept=x_intercept, color="red") +
                    geom_hline(yintercept=y_intercept, color="red") +
                    facet_wrap(as.formula(paste("~", facet_by))) +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(color=guide_legend(legend_title)) +
                    ggtitle(plot_title)
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(
                    print(
                        ggplot(data, aes_string(x=x_axis, y=y_axis, color=color_by)) +
                        geom_point(alpha=alpha) +
                        scale_colour_gradient(low="gray90", high="black") +
                        stat_smooth(method=lm) +
                        scale_x_log10() + 
                        scale_y_log10() + 
                        geom_vline(xintercept=x_intercept, color="red") +
                        geom_hline(yintercept=y_intercept, color="red") +
                        facet_wrap(as.formula(paste("~", facet_by))) +
                        xlab(x_label) +
                        ylab(y_label) +
                        guides(color=guide_legend(legend_title)) +
                        ggtitle(plot_title)
                    )
                )
                dev.off()
            }
            print(paste("Export log10 point plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export log10 point plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_vln_plot <- function(data, features, rootname, group_by=NULL, pt_size=NULL, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(
                print(
                    VlnPlot(
                        data,
                        features=features,
                        pt.size = pt_size,
                        group.by=group_by
                    )
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(
                    print(
                        VlnPlot(
                            data,
                            features=features,
                            pt.size = pt_size,
                            group.by=group_by
                        )
                    )
                )
                dev.off()
            }
            print(paste("Export violin plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export violin plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_feature_plot <- function(data, features, rootname, reduction, split_by=NULL, label=FALSE, order=FALSE, min_cutoff=NA, max_cutoff=NA, pt_size=NULL, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(
                print(
                    FeaturePlot(
                        data,
                        features=features,
                        pt.size=pt_size,
                        order=order,
                        min.cutoff=min_cutoff,
                        max.cutoff=max_cutoff,
                        reduction=reduction,
                        split.by=split_by,
                        label=label
                    )
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(
                    print(
                        FeaturePlot(
                            data,
                            features=features,
                            pt.size=pt_size,
                            order=order,
                            min.cutoff=min_cutoff,
                            max.cutoff=max_cutoff,
                            reduction=reduction,
                            split.by=split_by,
                            label=label
                        )
                    )
                )
                dev.off()
            }
            print(paste("Export feature plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export feature plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_variable_feature_plot <- function(data, rootname, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(VariableFeaturePlot(data)))
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(VariableFeaturePlot(data)))
                dev.off()
            }
            print(paste("Export variable feature plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export variable feature plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_dim_plot <- function(data, rootname, reduction, split_by=NULL, group_by=NULL, label=FALSE, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(
                print(
                    DimPlot(
                        data,
                        reduction=reduction,
                        split.by=split_by,
                        group.by=group_by,
                        label=label
                    )
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(
                    print(
                        DimPlot(
                            data,
                            reduction=reduction,
                            split.by=split_by,
                            group.by=group_by,
                            label=label
                        )
                    )
                )
                dev.off()
            }
            print(paste("Export Dim plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export Dim plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_pca_heatmap <- function(data, rootname, dims=NULL, nfeatures=30, cells=500, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(
                print(
                    DimHeatmap(
                        data,
                        dims=dims,
                        nfeatures=nfeatures,
                        reduction="pca",
                        cells=cells,
                        balanced=TRUE
                    )
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(
                    print(
                        DimHeatmap(
                            data,
                            dims=dims,
                            nfeatures=nfeatures,
                            reduction="pca",
                            cells=cells,
                            balanced=TRUE
                        )
                    )
                )
                dev.off()
            }
            print(paste("Export PCA heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export PCA heatmap to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_pca_loadings_plot <- function(data, rootname, dims=NULL, nfeatures=30, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(
                print(
                    VizDimLoadings(
                        data,
                        dims=dims,
                        nfeatures=nfeatures,
                        reduction="pca"
                    )
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(
                    print(
                        VizDimLoadings(
                            data,
                            dims=dims,
                            nfeatures=nfeatures,
                            reduction="pca"
                        )
                    )
                )
                dev.off()
            }
            print(paste("Export PCA loadings plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export PCA loadings plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_elbow_plot <- function(data, rootname, ndims=NULL, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(
                print(
                    ElbowPlot(
                        data,
                        ndims=ndims,
                        reduction="pca"
                    )
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(
                    print(
                        ElbowPlot(
                            data,
                            ndims=ndims,
                            reduction="pca"
                        )
                    )
                )
                dev.off()
            }
            print(paste("Export Elbow plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            dev.off()
            print(paste("Failed to export Elbow plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_all_clustering_plots <- function(seurat_data, suffix, args) {
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        export_dim_plot(
            data=seurat_data,
            reduction="umap",
            split_by="condition",
            group_by=paste("integrated_snn_res", current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "umap_plot_spl_by_cond_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="umap",
            group_by=paste("cluster_ext_type_res", current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "umap_ctype_pred_plot_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="umap",
            split_by="Phase",
            group_by=paste("integrated_snn_res", current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "umap_plot_spl_by_ph_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        Idents(seurat_data) <- paste("integrated_snn_res", current_resolution, sep=".")
        export_feature_plot(
            data=seurat_data,
            features=c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "mito_percentage", "log10_gene_per_log10_umi"),
            reduction="umap",
            label=TRUE,
            order=TRUE,
            min_cutoff="q10",
            rootname=paste(args$output, suffix, "umap_qc_metrics_plot_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        Idents(seurat_data) <- "new.ident"
    }
}


export_all_dimensionality_plots <- function(seurat_data, suffix, args) {
    export_elbow_plot(
        data=seurat_data,
        ndims=50,
        rootname=paste(args$output, suffix, "elbow_plot", sep="_"),
        pdf=args$pdf
    )
    export_dim_plot(
        data=seurat_data,
        reduction="pca",
        rootname=paste(args$output, suffix, "pca_plot", sep="_"),
        pdf=args$pdf
    )
    export_pca_heatmap(
        data=seurat_data,
        dims=1:50,
        rootname=paste(args$output, suffix, "pca_heatmap", sep="_"),
        height=4200,
        pdf=args$pdf
    )
    export_pca_loadings_plot(
        data=seurat_data,
        dims=1:50,
        rootname=paste(args$output, suffix, "pca_loadings_plot", sep="_"),
        height=4200,
        pdf=args$pdf
    )
    export_dim_plot(
        data=seurat_data,
        reduction="umap",
        split_by="new.ident",
        rootname=paste(args$output, suffix, "umap_plot_spl_by_ident", sep="_"),
        pdf=args$pdf
    )
}


export_all_qc_plots <- function(seurat_data, suffix, args){
    export_geom_bar_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "cell_count_plot", sep="_"),
        x_axis="new.ident",
        color_by="new.ident",
        x_label="Identity",
        y_label="Cells",
        legend_title="Identity",
        plot_title=paste("Number of cells per dataset (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_log10_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "umi_log10_density_plot", sep="_"),
        x_axis="nCount_RNA",
        color_by="new.ident",
        x_intercept=args$minumi,
        x_label="UMI's per cell",
        y_label="Density (log10)",
        legend_title="Identity",
        plot_title=paste("UMI density per cell (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_log10_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gene_log10_density_plot", sep="_"),
        x_axis="nFeature_RNA",
        color_by="new.ident",
        x_intercept=args$minfeatures,
        x_label="Genes per cell",
        y_label="Density (log10)",
        legend_title="Identity",
        plot_title=paste("Gene density per cell (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_log10_geom_point_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gene_umi_log10_correlation_plot", sep="_"),
        x_axis="nCount_RNA",
        y_axis="nFeature_RNA",
        color_by="mito_percentage",
        facet_by="new.ident",
        x_intercept=args$minumi,
        y_intercept=args$minfeatures,
        x_label="UMI's (log10)",
        y_label="Genes (log10)",
        legend_title="Mitochondrial %",
        plot_title=paste("Gene vs UMI correlation (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_log10_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "mito_perc_log10_density_plot", sep="_"),
        x_axis="mito_percentage",
        color_by="new.ident",
        x_intercept=args$maxmt,
        x_label="Mitochondrial gene percentage per cell",
        y_label="Density (log10)",
        legend_title="Identity",
        plot_title=paste("Mitochondrial gene percentage density per cell (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_log10_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "novelty_score_log10_density_plot", sep="_"),
        x_axis="log10_gene_per_log10_umi",
        color_by="new.ident",
        x_intercept=args$minnovelty,
        x_label="log10 Gene / log10 UMI per cell",
        y_label="Density (log10)",
        legend_title="Identity",
        plot_title=paste("Novelty score (log10Gene/log10UMI) density per cell (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_vln_plot(
        data=seurat_data,
        features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi"),
        rootname=paste(args$output, suffix, "qc_metrics_vln_plot", sep="_"),
        pt_size=0,
        pdf=args$pdf
    )
    export_vln_plot(
        data=seurat_data,
        features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi"),
        rootname=paste(args$output, suffix, "qc_metrics_vln_plot_gr_by_cond", sep="_"),
        group_by="condition",
        pt_size=0,
        pdf=args$pdf
    )
}


export_cellbrowser_data <- function(seurat_data, assay, matrix_slot, resolution, rootname){
    tryCatch(
        expr = {
            backup_assay <- DefaultAssay(seurat_data)
            DefaultAssay(seurat_data) <- assay
            ExportToCellbrowser(
                seurat_data,
                matrix.slot=matrix_slot,
                dir=rootname,
                cluster.field="Identity",
                meta.fields=c(
                    c("new.ident", "condition", "nCount_RNA", "nFeature_RNA", "log10_gene_per_log10_umi", "mito_percentage", "Phase", "S.Score", "G2M.Score"),
                    paste("integrated_snn_res", resolution, sep="."),
                    paste("cluster_ext_type_res", resolution, sep=".")
                ),
                meta.fields.names=c(
                    c("Identity", "Condition", "UMIs/cell", "Genes/cell", "Novelty score", "Mitochondrial %", "Cell Cycle Phase", "S phase score", "G2M phase score"),
                    paste("Clustering (", resolution, ")", sep=""),
                    paste("Cell Type Prediction (", resolution, ")", sep="")
                )
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


export_data <- function(data, location, row_names=FALSE, col_names=TRUE, quote=FALSE){
    tryCatch(
        expr = {
            write.table(
                data,
                file=location,
                sep=get_file_type(location),
                row.names=row_names,
                col.names=col_names,
                quote=quote
            )
            print(paste("Export data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data to", location, sep=" "))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description='Runs Seurat for comparative scRNA-seq analysis of across experimental conditions')
    # Import data from Cellranger Aggregate results
    parser$add_argument("--mex",           help="Path to the folder with not normalized aggregated feature-barcode matrices in MEX format", type="character", required="True")
    parser$add_argument("--identity",      help="Path to the aggregation CSV file to set the initial cell identity classes", type="character", required="True")
    parser$add_argument("--cellcycle",     help="Path to the TSV/CSV file with cell cycle data. First column - 'phase', second column 'gene_id'", type="character", required="True")
    parser$add_argument("--condition",     help="Path to the TSV/CSV file to define datasets conditions for grouping. First column - 'library_id' with values from the --identity file, second column 'condition'. Default: each dataset is assigned to its own biological condition", type="character")
    parser$add_argument("--classifier",    help="Path to the Garnett classifier rds file for cell type prediction", type="character", required="True")
    parser$add_argument("--barcodes",      help="Path to the headerless TSV/CSV file with selected barcodes (one per line) to prefilter input feature-barcode matrices. Default: use all cells", type="character")
    # Apply QC filters
    parser$add_argument("--mincells",      help="Include features detected in at least this many cells (applied to thoughout all datasets together). Default: 10", type="integer", default=10)
    parser$add_argument("--minfeatures",   help="Include cells where at least this many features are detected. Default: 250", type="integer", default=250)
    parser$add_argument("--minumi",        help="Include cells where at least this many UMI are detected. Default: 500", type="integer", default=500)
    parser$add_argument("--minnovelty",    help="Include cells with the novelty score not lower that this value (calculated as log10(genes)/log10(UMIs)). Default: 0.8", type="double", default=0.8)
    parser$add_argument("--maxmt",         help="Include cells with the mitochondrial contamination percentage not bigger that this value. Default: 5", type="double", default=5)
    parser$add_argument("--mitopattern",   help="Regex pattern to identify mitochondrial reads. Default: ^Mt-", type="character", default="^Mt-")
    # Integration, clustering, and cell types and marker genes identification parameters
    parser$add_argument("--regresscellcycle", help="Regress cell cycle as a confounding source of variation. Default: false", action="store_true")
    parser$add_argument("--regressmt",     help="Regress mitochondrial gene expression as a confounding source of variation. Default: false", action="store_true")
    parser$add_argument("--highvarcount",  help="Number of higly variable features to detect. Default: 3000", type="integer", default=3000)
    parser$add_argument("--ndim",          help="Number of principal components to use in clustering (1:50). Use Elbow plot to adjust this parameter. Default: 10", type="integer", default=10)
    parser$add_argument("--resolution",    help="Clustering resolution. Can be set as array. Default: 0.4 0.6 0.8 1.0 1.4", type="double", default=c(0.4, 0.6, 0.8, 1.0, 1.4), nargs="*")
    parser$add_argument("--logfc",         help="Log fold change threshold for conserved gene markers identification. Default: 0.25", type="double", default=0.25)
    parser$add_argument("--minpct",        help="Minimum fraction of cells where genes used for conserved gene markers identification should be detected in either of two tested clusters. Default: 0.1", type="double", default=0.1)
    parser$add_argument("--onlypos",       help="Return only positive markers when running conserved gene markers identification. Default: false", action="store_true")
    parser$add_argument("--testuse",       help="Set test type to use for putative and conserved gene marker identification. Default: wilcox", type="character", default="wilcox", choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"))
    parser$add_argument("--species",       help="Select species for gene name conversion when running cell type prediction", type="character", choices=names(SPECIES_DATA), required="True")
    # Export results
    parser$add_argument("--pdf",           help="Export plots in PDF. Default: false", action="store_true")
    parser$add_argument("--rds",           help="Save Seurat data to RDS file. Default: false", action="store_true")
    parser$add_argument("--output",        help="Output prefix. Default: ./seurat", type="character", default="./seurat")
    # Performance parameters
    parser$add_argument("--threads",       help="Threads. Default: 1", type="integer", default=1)
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}


args <- get_args()
print(paste("Setting parallelizations threads to", args$threads))
set_threads(args$threads)
print(paste("Loading cell identity data from Cellranger aggregation metadata file", args$identity))
cell_identity_data <- load_cell_identity_data (args$identity)
print(paste("Loading cell cycle data from", args$cellcycle))
cell_cycle_data <- load_cell_cycle_data(args$cellcycle)
print("Trying to load condition data")
condition_data <- load_condition_data(args$condition, cell_identity_data)
print(paste("Loading feature-barcode matrices from", args$mex))
seurat_data <- load_seurat_data(args$mex, args$mincells, cell_identity_data, condition_data)
print("Trying to load barcodes of interest to prefilter feature-barcode matrices by cells")
barcodes_data <- load_barcodes_data(args$barcodes, seurat_data)
print("Prefiltering feature-barcode matrices by cells of interest")
seurat_data <- apply_cell_filters(seurat_data, barcodes_data)
print(paste("Loading Garnett classifier from", args$classifier))
classifier <- readRDS(args$classifier)
print("Adding QC metrics to not filtered seurat data")
seurat_data <- add_qc_metrics(seurat_data, args$mitopattern)
export_all_qc_plots(seurat_data, "raw", args)
print("Applying QC filters to all datasets at once")
seurat_data <- apply_qc_filters(seurat_data, args$minfeatures, args$minumi, args$minnovelty, args$maxmt)
export_all_qc_plots(seurat_data, "filt", args)
print("Evaluating effects of cell cycle and mitochodrial gene expression")
explore_unwanted_variation(seurat_data, cell_cycle_data, args)
print("Running dataset integration")
seurat_data <- integrate_seurat_data(seurat_data, cell_cycle_data, args)
print("Performing PCA reduction of integrated data. Use all 50 principal components")
seurat_data <- RunPCA(seurat_data, npcs=50, verbose=FALSE)
print(paste("Performing UMAP reduction of integrated data using", args$ndim, "principal components"))
seurat_data <- RunUMAP(seurat_data, reduction="pca", dims=1:args$ndim, verbose=FALSE)
export_all_dimensionality_plots(seurat_data, "filt_int", args)
print(paste("Clustering integrated data using", args$ndim, "principal components"))
seurat_data <- FindNeighbors(seurat_data, reduction="pca", dims=1:args$ndim)
seurat_data <- FindClusters(seurat_data, resolution=args$resolution)
print("Assigning cell types for all clusters and all resolutions using only highly variable genes")
seurat_data <- assign_cell_types(seurat_data, classifier, args)
export_all_clustering_plots(seurat_data, "filt_int_cl", args)
export_rds(seurat_data, paste(args$output, "_data.rds", sep=""))
print("Identifying putative gene markers for all clusters and all resolutions")
all_putative_markers <- get_all_putative_markers(seurat_data, args)
export_data(all_putative_markers, paste(args$output, "_putative_gene_markers.tsv", sep=""))
print("Identifying conserved markers for all clusters and all resolutions irrespective of condition")
all_conserved_markers <- get_all_conserved_markers(seurat_data, args)
export_data(all_conserved_markers, paste(args$output, "_conserved_gene_markers.tsv", sep=""))
print("Exporting UCSC Cellbrowser data")
export_cellbrowser_data(seurat_data, assay="RNA", matrix_slot="data", resolution=args$resolution, rootname=paste(args$output, "_cellbrowser", sep=""))
