#!/usr/bin/env Rscript
options(warn=-1)
options(width=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(sva))
suppressMessages(library(cmapR))
suppressMessages(library(dplyr))
suppressMessages(library(future))
suppressMessages(library(ggpubr))
suppressMessages(library(hopach))
suppressMessages(library(Glimma))
suppressMessages(library(MAnorm2))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(argparse))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(htmlwidgets))
suppressMessages(library(BiocParallel))
suppressMessages(library(EnhancedVolcano))

theme_set(theme_classic())                       # set classic theme for all ggplot generated graphics

# To read:
# 1. https://cran.r-project.org/web/packages/MAnorm2/vignettes/MAnorm2_vignette.html

D40_COLORS <- c("#FB1C0D", "#0DE400", "#0D00FF", "#E8B4BD", "#FD00EA", "#0DD1FE", "#FF9B0D", "#0D601C", "#C50D69", "#CACA16", "#722A91", "#00DEBF", "#863B00", "#5D7C91", "#FD84D8", "#C100FB", "#8499FC", "#FD6658", "#83D87A", "#968549", "#DEB6FB", "#832E60", "#A8CAB0", "#FE8F95", "#FE1CBB", "#DF7CF8", "#FF0078", "#F9B781", "#4D493B", "#1C5198", "#7C32CE", "#EFBC16", "#7CD2DE", "#B30DA7", "#9FC0F6", "#7A940D", "#9B0000", "#946D9B", "#C8C2D9", "#94605A")

parallel <- function (args) {
    print(
        paste(
            "Setting parallelization to", args$cpus, "cores, and", args$memory,
            "GB of memory allowed to be shared between the processes"
        )
    )
    invisible(utils::capture.output(future::plan("multicore", workers=args$cpus)))
    invisible(utils::capture.output(future::plan()))
    invisible(utils::capture.output(data.table::setDTthreads(args$cpus)))
    options(future.globals.maxSize = args$memory * 1024^3)
    BiocParallel::register(BiocParallel::MulticoreParam(args$cpus, RNGseed=args$seed))
    set.seed(args$seed)
}

get_file_type <- function (filename){
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
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
            print(paste("Exporting data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data to", location, "due to", e))
        }
    )
}

export_gct <- function(counts_mat, row_metadata, col_metadata, location){
    tryCatch(
        expr = {
            row_metadata <- row_metadata %>% rownames_to_column("id") %>% mutate_at("id", as.vector)
            col_metadata <- col_metadata %>% rownames_to_column("id") %>% mutate_at("id", as.vector)
            gct_data <- new(
                "GCT",
                mat=counts_mat[row_metadata$id, col_metadata$id],       # to guarantee the order and number of row/columns
                rdesc=row_metadata,
                cdesc=col_metadata
            )
            write_gct(
                ds=gct_data,
                ofile=location,
                appenddim=FALSE
            )
            print(paste("Exporting GCT data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export GCT data to", location, sep=" "))
        }
    )
}

export_volcano_plot <- function(
        data, rootname,
        x_axis, y_axis,
        x_cutoff, y_cutoff,
        x_label, y_label,
        plot_title, plot_subtitle, caption,
        features=NA, pdf=FALSE, width=1200, height=800, resolution=100
){
    tryCatch(
        expr = {
            plot <- EnhancedVolcano(
                        data,
                        x=x_axis,
                        y=y_axis,
                        lab=rownames(data),
                        selectLab=features,
                        FCcutoff=x_cutoff,
                        pCutoff=y_cutoff,
                        xlab=x_label,
                        ylab=y_label,
                        title=plot_title,
                        subtitle=plot_subtitle,
                        caption=caption,
                        labSize=4,
                        labFace="bold",
                        labCol="red4",
                        colAlpha=0.6,
                        col=c("darkblue", "darkblue", "darkblue", "darkred"),
                        drawConnectors=TRUE,
                        widthConnectors=0.75
                    ) +
                    scale_y_log10() + annotation_logticks(sides="l", alpha=0.3) +
                    theme_classic() +
                    theme(legend.position="none", plot.subtitle=element_text(size=8, face="italic", color="gray30"))

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export volcano plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){})
            print(paste("Failed to export volcano plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

export_ma_plot <- function(
        data, rootname,
        fdr_cutoff, y_cutoff,
        x_label, y_label,
        plot_title, plot_subtitle, caption,
        pdf=FALSE, width=800, height=800, resolution=100
){
    tryCatch(
        expr = {
            plot <- ggmaplot(data,
                        fdr=fdr_cutoff,
                        fc=y_cutoff,
                        top=0,
                        size=2,
                        alpha=0.6,
                        main=plot_title,
                        subtitle=plot_subtitle,
                        caption=caption,
                        xlab=x_label,
                        ylab=y_label,
                        palette=c("darkred", "darkred", "darkblue")
                    ) +
                    theme_classic() +
                    theme(legend.position="none", plot.subtitle=element_text(size=8, face="italic", color="gray30"))

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (!is.null(pdf) && pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Exporting MA-plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){})
            print(
                paste0(
                    "Failed to export MA-plot to ", rootname,
                    ".(png/pdf) with error - ", e
                )
            )
        }
    )
}

get_pca_data <- function(counts_data){
    tryCatch(
        expr = {
            print("Computing PCA for counts data")
            counts_data <- counts_data %>%
                           filter_all(any_vars(. != 0))   # remove rows with only zeros, otherwise prcomp may fail
            pca_raw <- prcomp(
                t(counts_data),
                center=TRUE,
                scale.=TRUE
            )
            pca_scores <- as.data.frame(pca_raw$x) %>%
                          rownames_to_column(var="group")
            pca_variance <- round(pca_raw$sdev / sum(pca_raw$sdev) * 100, 2)
            return (list(scores=pca_scores, variance=pca_variance))
        },
        error = function(e){
            print(paste("Failed to compute PCA for counts data due to", e))
        }
    )
}

export_pca_plot <- function(
        data, pcs, rootname, plot_title, legend_title,
        plot_subtitle=NULL, color_by="label",
        label_size=5, pt_size=8, pt_shape=19, alpha=0.75,
        palette_colors=D40_COLORS,
        pdf=FALSE, width=1200, height=800, resolution=100
){
    tryCatch(
        expr = {
            x_score_column <- paste0("PC", pcs[1])
            y_score_column <- paste0("PC", pcs[2])
            x_variance <- data$variance[pcs[1]]
            y_variance <- data$variance[pcs[2]]
            plot <- ggplot(
                        data$scores,
                        aes_string(x=x_score_column, y=y_score_column, color=color_by)
                    ) +
                    geom_point(size=pt_size, shape=pt_shape, alpha=alpha) +
                    xlab(paste0(x_score_column, ": ", x_variance, "% variance")) +
                    ylab(paste0(y_score_column, ": ", y_variance, "% variance")) + 
                    geom_label_repel(
                        aes_string(label=color_by),
                        size=label_size,
                        point.padding=0.5,
                        box.padding=0.5,
                        check_overlap=TRUE,
                        show.legend=FALSE
                    ) +
                    ggplot2::ggtitle(plot_title, subtitle=plot_subtitle) +
                    guides(color=guide_legend(legend_title)) +
                    scale_color_manual(values=palette_colors) +

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (!is.null(pdf) && pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Exporting PCA plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){})
            print(paste("Failed to export PCA plot to ", rootname, ".(png/pdf) with error - ", e, sep=""))
        }
    )
}

export_mds_html_plot <- function(data, groups, location){
    tryCatch(
        expr = {
            saveWidget(
                glimmaMDS(
                    x=data,
                    groups=groups,
                    labels=rownames(groups)
                ),
                file=location
            )
        },
        error = function(e){
            print(paste0("Failed to export MDS plot to ", location, " with error - ", e))
        }
    )
}

export_plots <- function(db_sites, metadata, args){
    export_volcano_plot(
        data=db_sites,
        x_axis="log2FoldChange",
        y_axis="padj",
        x_cutoff=0,
        y_cutoff=args$padj,
        x_label=bquote(~Log[2]~"fold change"),
        y_label=bquote("-"~Log[10]~"padj"),
        plot_title="Volcano plot for differentially bound sites",
        plot_subtitle=paste0("Differentially bound sites with padj", "<=", args$padj),
        caption=paste(nrow(db_sites), "features"),
        rootname=paste(args$output, "diff_vlcn", sep="_"),
        pdf=args$pdf
    )

    export_ma_plot(
        data=db_sites,
        fdr_cutoff=args$padj,
        y_cutoff=0,
        x_label=bquote(~Log[2]~"base mean"),
        y_label=bquote(~Log[2]~"fold change"),
        plot_title="MA-plot for differentially bound sites",
        plot_subtitle=paste0("Differentially bound sites with padj", "<=", args$padj),
        caption=paste(nrow(data), "features"),
        rootname=paste(args$output, "diff_ma", sep="_"),
        pdf=args$pdf
    )

    pca_data <- get_pca_data(db_sites[, rownames(metadata)])               # adds 'group' column to identify the datasets

    export_pca_plot(
        data=pca_data,
        pcs=c(1, 2),
        plot_title="Read counts PCA",
        plot_subtitle="PC1/PC2",
        legend_title="Dataset",
        color_by="group",
        width=ifelse(length(rownames(metadata)) > 13, 1600, 1200),       # need to make it wider if we have a lot of samples
        rootname=paste(args$output, "pca_1_2", sep="_"),
        pdf=args$pdf
    )

    export_pca_plot(
        data=pca_data,
        pcs=c(2, 3),
        plot_title="Read counts PCA",
        plot_subtitle="PC2/PC3",
        legend_title="Dataset",
        color_by="group",
        width=ifelse(length(rownames(metadata)) > 13, 1600, 1200),       # need to make it wider if we have a lot of samples
        rootname=paste(args$output, "pca_2_3", sep="_"),
        pdf=args$pdf
    )

    export_mds_html_plot(
        data=as.matrix(db_sites[, rownames(metadata)]),
        groups=metadata,
        location=paste(args$output, "mds_plot.html", sep="_")
    )

}

assert_args <- function(args){
    if ((length(args$read1) != length(args$peak1)) || (length(args$read1) != length(args$name1))){
        print(
            paste(
                "Number of the values provided for the --read1,",
                "--peak1, and --name1 parameters doesn't match"
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
    if ((length(args$read2) != length(args$peak2)) || (length(args$read2) != length(args$name2))){
        print(
            paste(
                "Number of the values provided for the --read2,",
                "--peak2, and --name2 parameters doesn't match"
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
    if(!is.null(args$summit1) && (length(args$summit1) != length(args$read1))){
        print(
            paste(
                "Number of the values provided for the --summit1",
                "and --read1 parameters doesn't match"
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
    if(!is.null(args$summit2) && (length(args$summit2) != length(args$read2))){
        print(
            paste(
                "Number of the values provided for the --summit2",
                "and --read2 parameters doesn't match"
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
    if(xor(is.null(args$summit1), is.null(args$summit2))){
        print(
            paste(
                "Either both or neither of the --summit1 and",
                "--summit2 parameters should be provided"
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
    if(any(duplicated(c(args$name1, args$name2)))){
        print(
            paste(
                "--name1 and/or --name2 parameters",
                "include repeated values"
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
    args$minoverlap1 <- ifelse(
        as.integer(args$minoverlap) == args$minoverlap,
        min(args$minoverlap, length(args$name1)),
        max(round(args$minoverlap * length(args$name1)), 1)
    )
    args$minoverlap2 <- ifelse(
        as.integer(args$minoverlap) == args$minoverlap,
        min(args$minoverlap, length(args$name2)),
        max(round(args$minoverlap * length(args$name2)), 1)
    )
    return (args)
}

export_corr_heatmap <- function(
        data, rootname, plot_title, x_label, y_label,
        plot_subtitle=NULL, data_labels=NULL,
        font_color="white", font_size=4, fontface="bold",
        gradient_colors=c("darkblue", "lightgrey", "darkred"),
        pdf=FALSE, width=800, height=800, resolution=100
){
    base::tryCatch(
        expr = {
            if (!is.null(data_labels)){
                colnames(data) <- data_labels
                rownames(data) <- data_labels
            }
            data_long <- melt(data)
            plot <- ggplot(
                        data_long, aes(Var1, Var2, fill=value)
                    ) +
                    geom_tile() +
                    geom_text(
                        aes(label=round(value, 2)),
                        color=font_color,
                        size=font_size,
                        fontface=fontface
                    ) +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(
                        x=guide_axis(angle=45),
                        y=guide_axis(angle=45)
                    ) +
                    ggtitle(plot_title, subtitle=plot_subtitle) +
                    scale_fill_gradient2(
                        low=gradient_colors[1], mid=gradient_colors[2], high=gradient_colors[3],
                        midpoint=0, limit=c(-1, 1),
                        name="Correlation"
                    ) +
                    coord_fixed(ratio=1)

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(plot))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(plot))
                grDevices::dev.off()
            }

            base::print(base::paste("Exporting samples correlation heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(
                base::paste0(
                    "Failed to export samples correlation heatmap to ",
                    rootname, ".(png/pdf) with error - ", e
                )
            )
        }
    )
}

export_overlap_heatmap <- function(
        data, rootname, plot_title, x_label, y_label,
        plot_subtitle=NULL, data_labels=NULL,
        font_color="white", font_size=4, fontface="bold",
        overlap_colors=c("lightgrey", "darkblue"),
        jaccard_colors=c("lightgrey", "darkblue"),
        combine_guides="keep", pdf=FALSE, width=1200, height=800, resolution=100
){
    base::tryCatch(
        expr = {

            jaccard_index <- function(x, y) {
                intersection <- sum(x & y)
                union <- sum(x | y)
                return (intersection/union)
            }

            overlap_matrix <- t(as.matrix(data)) %*% as.matrix(data)
            if (!is.null(data_labels)){
                colnames(overlap_matrix) <- data_labels
                rownames(overlap_matrix) <- data_labels
            }
            overlap_matrix_long <- melt(overlap_matrix)

            jaccard_matrix <- matrix(0, nrow=ncol(data), ncol=ncol(data))
            for (i in 1:ncol(data)) {
                for (j in 1:ncol(data)) {
                    jaccard_matrix[i, j] <- jaccard_index(data[, i], data[, j])
                }
            }
            if (!is.null(data_labels)){
                colnames(jaccard_matrix) <- data_labels
                rownames(jaccard_matrix) <- data_labels
            }
            jaccard_matrix_long <- melt(jaccard_matrix)

            plot_overlap <- ggplot(
                                overlap_matrix_long, aes(Var1, Var2, fill=value)
                            ) +
                            geom_tile() +
                            geom_text(
                                aes(label=scales::comma(value)),
                                color=font_color,
                                size=font_size,
                                fontface=fontface
                            ) +
                            xlab(x_label) +
                            ylab(y_label) +
                            guides(
                                x=guide_axis(angle=45),
                                y=guide_axis(angle=45)
                            ) +
                            ggtitle("Overlaps") +
                            scale_fill_gradient(trans="log10", low=overlap_colors[1], high=overlap_colors[2], name="Overlap\nCounts") +
                            coord_fixed(ratio=1)

            plot_jaccard <- ggplot(
                                jaccard_matrix_long, aes(Var1, Var2, fill=value)
                            ) +
                            geom_tile() +
                            geom_text(
                                aes(label=round(value, 2)),
                                color=font_color,
                                size=font_size,
                                fontface=fontface
                            ) +
                            xlab(x_label) +
                            ylab(y_label) +
                            guides(
                                x=guide_axis(angle=45),
                                y=guide_axis(angle=45)
                            ) +
                            ggtitle("Jaccard") +
                            scale_fill_gradient(low=jaccard_colors[1], high=jaccard_colors[2], limit=c(0, 1), name="Jaccard\nIndex") +
                            coord_fixed(ratio=1)

            combined_plots <- wrap_plots(plot_overlap + plot_jaccard, guides=combine_guides) +
                              plot_annotation(
                                  title=plot_title,
                                  subtitle=plot_subtitle
                              )

            grDevices::png(filename=base::paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            base::suppressMessages(base::print(combined_plots))
            grDevices::dev.off()

            if (!is.null(pdf) && pdf) {
                grDevices::pdf(file=base::paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                base::suppressMessages(base::print(combined_plots))
                grDevices::dev.off()
            }

            base::print(base::paste("Exporting samples overlap heatmap to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            base::tryCatch(expr={grDevices::dev.off()}, error=function(e){})
            base::print(
                base::paste0(
                    "Failed to export samples overlap heatmap to ",
                    rootname, ".(png/pdf) with error - ", e
                )
            )
        }
    )
}

get_args <- function(){
    parser <- ArgumentParser(
        description="MAnorm2 for Normalizing and Comparing ChIP-seq Samples"
    )
    parser$add_argument(
        "--read1",
        help=paste(
            "Coordinate sorted and indexed BAM files with",
            "aligned reads from the samples that belong",
            "to the first biological condition."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--read2",
        help=paste(
            "Coordinate sorted and indexed BAM files with",
            "aligned reads from the samples that belong",
            "to the second biological condition."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--peak1",
        help=paste(
            "Narrow or broad peak files with the peaks",
            "called from the samples that belong to the",
            "first biological condition."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--peak2",
        help=paste(
            "Narrow or broad peak files with the peaks",
            "called from the samples that belong to the",
            "second biological condition."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--summit1",
        help=paste(
            "BED files with the summits of the peaks called",
            "from the samples that belong to the first",
            "biological condition. If not provided, the peak",
            "center is taken as the summit."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--summit2",
        help=paste(
            "BED files with the summits of the peaks called",
            "from the samples that belong to the second",
            "biological condition. If not provided, the peak",
            "center is taken as the summit."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--name1",
        help=paste(
            "Names of the samples that belong to the first",
            "biological condition. All values from the --name1",
            "and --name2 parameters should be unique."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--name2",
        help=paste(
            "Names of the samples that belong to the second",
            "biological condition. All values from the --name1",
            "and --name2 parameters should be unique."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--condition1",
        help=paste(
            "Name for the first biological condition. The",
            "direction of comparison is always --condition2",
            "vs --condition1. Default: control."
        ),
        type="character", default="control"
    )
    parser$add_argument(
        "--condition2",
        help=paste(
            "Name for the second biological condition. The",
            "direction of comparison is always --condition2",
            "vs --condition1. Default: treatment."
        ),
        type="character", default="treatment"
    )
    parser$add_argument(
        "--minoverlap",
        help=paste(
            "Filtering threshold to keep only those reference",
            "genomic bins that are present in at least this",
            "many samples within the biological condition. If",
            "this threshold has a value between zero and one,",
            "only those peaks will be included that are present in",
            "at least this fraction of datasets. Default: 1"
        ),
        type="double", default=1
    )
    parser$add_argument(
        "--padj",
        help=paste(
            "Filtering threshold to report only differentially",
            "bound sites with adjusted P-value less than or equal",
            "to the provided value. Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--batch",
        help=paste(
            "Optional headerless TSV/CSV file for batch effect",
            "correction. First column should include values from",
            "the --name1 and --name2 parameters in any order.",
            "The second column should include group names to",
            "which each sample should be assigned. Default: do",
            "not apply batch correction."
        ),
        type="character"
    )
    parser$add_argument(
        "--maxpeaks",
        help=paste(
            "The maximum number of the most significant peaks",
            "to select from each peak file when constructing",
            "reference genomic bins. The top significant peaks",
            "are selected based on the score column which is",
            "calculated by MACS2 either as int(-10*log10qvalue)",
            "or as int(-10*log10qvalue), depending on the cutoff",
            "used for peak calling. Default: keep all peaks."
        ),
        type="integer"
    )
    parser$add_argument(
        "--minpeakgap",
        help=paste(
            "Peaks remained after optional filtering by --maxpeaks",
            "parameter will be merged if the distance between them",
            "is smaller than the provided value. Merging is first",
            "applied per sample and then to all peaks together",
            "before splitting them into the reference genomic bins",
            "of size --binsize. Default: 150"
        ),
        type="integer", default=150
    )
    parser$add_argument(
        "--binsize",
        help=paste(
            "The size of non-overlapping reference genomic bins used",
            "for generating a table of read per peak counts. 2000 bp",
            "is recommended for sharp histone marks like H3K4me3 and",
            "H3K27ac, and 1000 bp for TFs or DNase-seq. Default: 2000"
        ),
        type="integer", default=2000
    )
    parser$add_argument(
        "--fixbinsize",
        help=paste(
            "Force all reference genomic bins be exaclty the",
            "same size as provided in the --binsize parameter.",
            "Default: when a merged peak is split into the",
            "reference genomic bins, their sizes do not exceed",
            "the boundaries of the original peak."
        ),
        action="store_true"
    )
    parser$add_argument(
        "--blacklist",
        help=paste(
            "BED file with the genomic blacklist regions. Any",
            "reference genomic bin overlapping a blacklist region",
            "will be removed from the analysis. Default: include",
            "all reference genomic bins."
        ),
        type="character"
    )
    parser$add_argument(
        "--dedup",
        help=paste(
            "Remove duplicated reads identified by their coordinates.",
            "The location of a single-end read is determined by its",
            "strand and 5' coordinate. For a paired-end read, the",
            "DNA fragment position is used instead. Default: include",
            "all reads."
        ),
        action="store_true"
    )
    parser$add_argument(
        "--shiftsize",
        help=paste(
            "Shift the positions of the 5' ends of a single-end",
            "reads downstream on the selected value. Use the resulting",
            "points for counting reads in the reference genomic bins.",
            "Set as half of the DNA fragment size. Ignored if --paired",
            "parameter is provided. Default 100"
        ),
        type="integer", default=100
    )
    parser$add_argument(
        "--paired",
        help=paste(
            "Consider all reads as paired-end. When counting reads",
            "in the reference genomic bins, use the middle point of",
            "each DNA fragment. --shiftsize parameters is ignored.",
            "Default: treat all reads as single-end."
        ),
        action="store_true"
    )
    parser$add_argument(
        "--exclude",
        help=paste(
            "Define the chromosomes to be exluded from the analysis.",
            "Default: include all chromosomes"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--norm",
        help=paste(
            "Normalization method applied to the raw read counts.",
            "pseudo-reference - normalize each sample to the pseudo",
            "reference that includes the average intesities from all",
            "samples. A reference genomic bin is occupied by the",
            "pseudo reference if it was occupied by at least one",
            "sample that the reference was constructed from. Each",
            "sample is MA-normalized to the pseudo reference using",
            "the common genomic bins between the reference and a",
            "sample. baseline - normalize each sample to the one",
            "whose log2 size factor is closest to 0. hierarchical -",
            "similar to the baseline but first all samples are",
            "normalized within the biological conditions, than two",
            "biological conditions are normalized between each other.",
            "Default: pseudo-reference"
        ),
        type="character",
        default="pseudo-reference",
        choices=c("pseudo-reference", "baseline", "hierarchical")
    )
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./manorm",
        type="character", default="./manorm"
    )
    parser$add_argument(
        "--cpus",
        help="Number of cores/cpus to use. Default: 1",
        type="integer", default=1
    )
    parser$add_argument(
        "--memory",
        help=paste(
            "Maximum memory in GB allowed to be shared between",
            "the workers when using multiple --cpus.",
            "Default: 32"
        ),
        type="integer", default=32
    )
    parser$add_argument(
        "--tmpdir",
        help=paste(
            "Directory to keep temporary files. Default: either",
            "/tmp or defined by environment variables TMPDIR,",
            "TMP, TEMP."
        ),
        type="character", default=tempdir()
    )
    parser$add_argument(
        "--seed",
        help=paste(
            "Seed number for random values. Default: 42"
        ),
        type="integer", default=42
    )
    args <- assert_args(parser$parse_args(commandArgs(trailingOnly = TRUE)))
    print(args)
    return (args)
}

args <- get_args()
parallel(args)
dir.create(args$tmpdir, recursive=TRUE)

metadata <- data.frame(
    bams=c(args$read1, args$read2),
    peaks=c(args$peak1, args$peak2),
    conditions=c(
        rep(args$condition1, length(args$name1)),
        rep(args$condition2, length(args$name2))
    ),
    check.names=FALSE
)
if(!is.null(args$summit1) && !is.null(args$summit2)){                      # we already checked in the assert_args that either both
    metadata$summits <- c(args$summit1, args$summit2)                      # or neither of the summit parameters are provided
}
rownames(metadata) <- c(args$name1, args$name2)
metadata$read_cnt <- paste(rownames(metadata), "read_cnt", sep=".")
metadata$occupancy <- paste(rownames(metadata), "occupancy", sep=".")
print(metadata)

all_beds <- c()
for (i in 1:nrow(metadata)){
    current_bam <- metadata$bams[i]
    current_sam <- file.path(args$tmpdir, paste0("sample_", i, ".sam"))
    current_bed <- file.path(args$tmpdir, paste0("sample_", i, ".bed"))

    print(paste("Attempting to convert", current_bam, "to", current_sam))
    exit_code <- sys::exec_wait(
        cmd="samtools",
        args=c(
            "view", "--with-header",
            "--output", current_sam,
            "--threads", args$cpus,
            current_bam
        )
    )
    if (exit_code != 0){
        print(paste("Failed to convert", current_bam, "to", current_sam))
        unlink(args$tmpdir, recursive=TRUE)
        quit(save="no", status=1, runLast=FALSE)
    }

    print(paste("Attempting to convert", current_sam, "to", current_bed))
    exit_code <- sys::exec_wait(
        cmd="sam2bed",
        args=c(
            "-i", current_sam,
            "-o", current_bed
        )
    )
    if (exit_code != 0){
        print(paste("Failed to convert", current_sam, "to", current_bed))
        unlink(args$tmpdir, recursive=TRUE)
        quit(save="no", status=1, runLast=FALSE)
    }
    file.remove(current_sam)
    all_beds <- c(all_beds, current_bed)
}
metadata$beds <- all_beds
print(metadata)

print("Counting reads in the reference genomic bins")
exit_code <- sys::exec_wait(
    cmd="profile_bins",
    args=na.omit(
        c(
            paste0("--reads=", paste(metadata$beds, collapse=",")),
            paste0("--peaks=", paste(metadata$peaks, collapse=",")),
            paste0("--labs=", paste(rownames(metadata), collapse=",")),
            ifelse(
                "summits" %in% colnames(metadata),
                paste0("--summits=", paste(metadata$summits, collapse=",")),
                NA
            ),
            ifelse(
                !is.null(args$maxpeaks),
                paste0("--keep-peaks=", args$maxpeaks),
                NA
            ),
            paste0("--min-peak-gap=", args$minpeakgap),
            ifelse(
                args$dedup,
                "--keep-dup=1",
                "--keep-dup=all"
            ),
            ifelse(args$paired, "--paired", NA),
            ifelse(args$fixbinsize, "--fix-bin-size", NA),
            paste0("--shiftsize=", args$shiftsize),
            paste0("--typical-bin-size=", args$binsize),
            ifelse(
                !is.null(args$blacklist),
                paste0("--filter=", args$blacklist),
                NA
            ),
            "-n", "peak"                                                             # this defines the prefix of the output file
        )
    )
)
if (exit_code != 0){
    print(paste("Failed to count reads in the reference genomic bins"))
    unlink(args$tmpdir, recursive=TRUE)
    quit(save="no", status=1, runLast=FALSE)
}

print("Loadind read counts data")
counts_data <- read.table(                                                           # region is "occupied" if its middle point is covered by a peak from a sample
    "peak_profile_bins.xls",
    sep="\t",
    header=TRUE,
    check.names=FALSE,
    stringsAsFactors=FALSE
)
if(!is.null(args$exclude)){
    print(
        paste(
            "Excluding", paste(args$exclude, collapse=", "),
            "chromosomes from the read counts matrix"
        )
    )
    print(paste("Before filtering", nrow(counts_data)))
    counts_data <- counts_data[!(counts_data$chrom %in% args$exclude),]
    print(paste("After filtering", nrow(counts_data)))
}
print(head(counts_data))

export_corr_heatmap(
    data=cor(counts_data[, metadata$read_cnt])[metadata$read_cnt, metadata$read_cnt],
    data_labels=rownames(metadata),
    x_label="Sample",
    y_label="Sample",
    plot_title="Read counts correlation between the samples",
    plot_subtitle="On the basis of the raw read counts\nwithin the reference genomic bins",
    rootname=paste(args$output, "smpl_corr_raw", sep="_"),
    pdf=args$pdf
)

if (!is.null(args$batch)){                                                           # need to correct batch effect
    print(
        paste(
            "Attempting to load batch metadata from", args$batch
        )
    )
    batch_data <- read.table(
        args$batch,
        sep=get_file_type(args$batch),
        row.names=1,
        col.names=c("sample", "batch"),                                              # "sample" columns is not used, but it fails without it 
        header=FALSE,
        stringsAsFactors=FALSE
    )
    if (!all(is.element(rownames(metadata), rownames(batch_data)))){                 # check if all our samples are present in the batch data file
        print(
            paste(
                "Batch metadata file has missing sample names. Exiting"
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
    metadata$batches <- batch_data[rownames(metadata), "batch"]                      # to make sure our batches are assigned in the right order
    print(metadata)

    print("Selecting all columns with the raw read counts")
    raw_data = counts_data[metadata$read_cnt]                                        # the order of columns will be the same as the rows in metadata

    print("Removing the batch effect using ComBat-Seq")
    corrected_data <- ComBat_seq(
        as.matrix(raw_data),
        batch=metadata$batches,                                                      # batch order will correspond to the columns in raw_data, because we selected them in the specific order
        group=metadata$conditions
    )
    counts_data[, colnames(corrected_data)] <- corrected_data                        # now all read_cnt columns have batch corrected data

    export_corr_heatmap(
        data=cor(counts_data[, metadata$read_cnt])[metadata$read_cnt, metadata$read_cnt],
        data_labels=rownames(metadata),
        x_label="Sample",
        y_label="Sample",
        plot_title="Read counts correlation between the samples",
        plot_subtitle="On the basis of the batch corrected raw read counts\nwithin the reference genomic bins",
        rootname=paste(args$output, "smpl_corr_crtd", sep="_"),
        pdf=args$pdf
    )
}
print(head(counts_data))

print(
    paste(
        "Normalizing read counts data using", args$norm, "mode"
    )
)
if (args$norm != "hierarchical"){
    counts_data <- MAnorm2::normalize(                                                # Constructs a pseudo ChIP-seq profile by “averaging” the intesities from all samples.
        counts_data,                                                                  # A reference genomic bin is occupied by the pseudo ChIP-seq sample if it was occupied
        count=metadata$read_cnt,                                                      # by at least one sample that it was constructed from. Then each sample is MA-normalized
        occupancy=metadata$occupancy,                                                 # to this pseudo reference using the common genomic bins between the reference and a sample.
        baseline=if (args$norm == "pseudo-reference") "pseudo-reference" else NULL
    )
} else {
    counts_data <- MAnorm2::normalize(
        MAnorm2::normalize(
            counts_data,
            count=metadata$read_cnt[metadata$conditions == args$condition1],
            occupancy=metadata$occupancy[metadata$conditions == args$condition1]
        ),
        count=metadata$read_cnt[metadata$conditions == args$condition2],
        occupancy=metadata$occupancy[metadata$conditions == args$condition2]
    )
}

print(
    paste0(
        "Adding datasets to the first (min ", args$minoverlap1,
        " overlap) and second (min ", args$minoverlap2,
        " overlap) biological conditions"
    )
)
bio_conditions <- list(                                                        # first, second, and common do not influence the columns names in the results
    first=MAnorm2::bioCond(
        norm.signal=counts_data[, metadata$read_cnt[metadata$conditions == args$condition1]],
        occupancy=counts_data[, metadata$occupancy[metadata$conditions == args$condition1]],
        name=args$condition1,                                                  # influences the name of the column in the results
        occupy.num=args$minoverlap1,
        meta.info=counts_data[, c("chrom", "start", "end")]                    # to have reference genomic bins coordinates embedded (not used by MAnorm2)
    ),
    second=MAnorm2::bioCond(
        norm.signal=counts_data[, metadata$read_cnt[metadata$conditions == args$condition2]],
        occupancy=counts_data[, metadata$occupancy[metadata$conditions == args$condition2]],
        name=args$condition2,                                                  # influences the name of the column in the results
        occupy.num=args$minoverlap2
    )
)
if (length(args$name1) == 1 && length(args$name2) == 1){                       # we have only two samples to compare, so need to add common condition
    print("Adding common biological condition")
    bio_conditions[["common"]] <- MAnorm2::bioCond(
        norm.signal=counts_data[metadata$read_cnt],
        occupancy=counts_data[metadata$occupancy],
        occupy.num=2,                                                          # should be always 2 because we have only 2 samples
        name="common"                                                          # influences the name of the column in the results
    )
}

if (args$norm == "hierarchical"){
    print("Performing between-group normalization")
    bio_conditions <- MAnorm2::normBioCond(bio_conditions)
    attr(counts_data, "size.factor") <- attr(bio_conditions, "size.factor")
    attr(counts_data, "baseline") <- attr(bio_conditions, "baseline")
    attr(counts_data, "norm.coef") <- attr(bio_conditions, "norm.coef")
    attr(counts_data, "MA.cor") <- attr(bio_conditions, "MA.cor")
}

counts_data[, colnames(bio_conditions$first$norm.signal)] <- bio_conditions$first$norm.signal
counts_data[, colnames(bio_conditions$second$norm.signal)] <- bio_conditions$second$norm.signal

print(head(counts_data))
print(attr(counts_data, "norm.coef"))
print(attr(counts_data, "MA.cor"))
print(attr(counts_data, "baseline"))

export_corr_heatmap(
    data=cor(counts_data[, metadata$read_cnt])[metadata$read_cnt, metadata$read_cnt],
    data_labels=rownames(metadata),
    x_label="Sample",
    y_label="Sample",
    plot_title="Read counts correlation between the samples",
    plot_subtitle=paste0(
        "On the basis of the ",
        ifelse(
            !is.null(args$batch),
            "batch corrected ",
            ""
        ),
        "normalized read counts\nwithin the reference genomic bins"
    ),
    rootname=paste(args$output, "smpl_corr_norm", sep="_"),
    pdf=args$pdf
)

export_overlap_heatmap(
    data=counts_data[, metadata$occupancy],                                          # occupancy has the same order as rownames, so data_labels will correspond to data
    data_labels=rownames(metadata),
    x_label="Sample",
    y_label="Sample",
    plot_title="Peaks overlap between the samples",
    plot_subtitle="On the basis of the occupied by each sample reference genomic bins",
    rootname=paste(args$output, "smpl_vrlp", sep="_"),
    pdf=args$pdf
)

export_overlap_heatmap(
    data=data.frame(
        first=as.numeric(bio_conditions$first$occupancy),                      # better to convert to numeric
        second=as.numeric(bio_conditions$second$occupancy),                    # better to convert to numeric
        check.names=FALSE
    ),
    data_labels=c(args$condition1, args$condition2),                           # data_labels order corresponds to the columns order in data input
    x_label="Condition",
    y_label="Condition",
    plot_title="Peaks overlap between the biological conditions",
    plot_subtitle="On the basis of the occupied by each biological condition reference genomic bins",
    rootname=paste(args$output, "cnd_vrlp", sep="_"),
    width=1200, height=600,
    pdf=args$pdf
)

if (args$norm == "hierarchical"){
    export_corr_heatmap(
        data=attr(counts_data, "MA.cor")[c(args$condition1, args$condition2), c(args$condition1, args$condition2)],
        data_labels=c(args$condition1, args$condition2),
        x_label="Condition",
        y_label="Condition",
        plot_title="Correlation between M and A values across\nthe common peak regions of each pair of biological conditions",
        plot_subtitle=paste("Normalized using", args$norm, "mode"),
        rootname=paste(args$output, "ma_corr", sep="_"),
        pdf=args$pdf
    )
} else {
    export_corr_heatmap(
        data=attr(counts_data, "MA.cor")[metadata$read_cnt, metadata$read_cnt],
        data_labels=rownames(metadata),
        x_label="Sample",
        y_label="Sample",
        plot_title="Correlation between M and A values across\nthe common peak regions of each pair of samples",
        plot_subtitle=paste("Normalized using", args$norm, "mode"),
        rootname=paste(args$output, "ma_corr", sep="_"),
        pdf=args$pdf
    )
}

print("Fitting mean-variance curve based on the common bins")
bio_conditions <- MAnorm2::fitMeanVarCurve(                                    # will add fit.info field to the bioCond objects
    bio_conditions,                                                            # at list one condition should have replicates 
    method="parametric",
    occupy.only=TRUE,
    init.coef=c(0.1, 10)                                                       # per manual it is expected to suit most practical datasets
)
print(bio_conditions)

print("Running differential test")
db_sites <- MAnorm2::diffTest(                                                       # data frame without coordinates, Mval is calculated as Y/X, not filtered
                x=bio_conditions$first,
                y=bio_conditions$second,
            ) %>%                                                                    # first.mean, second.mean, Mval, Mval.se, Mval.t, pval, padj
            tibble::rownames_to_column(var="rowname") %>%                            # we need it to join with coordinates
            dplyr::left_join(
                bio_conditions$first$meta.info %>%                                  # we take the coordinates from the meta.info of the second bio condition
                tibble::rownames_to_column(var="rowname"),
                by="rowname"
            ) %>%
            stats::na.omit() %>%                                                     # we shouldn't have any NAs, but filter just in case 
            dplyr::rename(
                "pvalue"="pval",
                "log2FoldChange"="Mval",
                "lfcSE"="Mval.se",                                                   # to have the name similar to what DESeq2 reports
                "padj"="padj",
                "chr"="chrom"
            ) %>%
            dplyr::mutate(
                "baseMean"=(.[[paste0(args$condition1, ".mean")]] + .[[paste0(args$condition2, ".mean")]]) / 2
            ) %>%
            dplyr::select(                                                           # to have a proper columns order
                c(
                    "chr", "start", "end",
                    "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj",
                    paste0(args$condition1, ".mean"),
                    paste0(args$condition2, ".mean")
                )
            ) %>% 
            dplyr::left_join(
                counts_data %>%
                    dplyr::select(-metadata$occupancy) %>%                           # keeping the occupancy columns might be misleading
                    dplyr::rename("chr"="chrom") %>%
                    dplyr::rename(!!!setNames(metadata$read_cnt, rownames(metadata))),
                by=c("chr", "start", "end")
            ) %>%
            dplyr::mutate(
                "feature"=paste(.$chr, paste(.$start, .$end, sep="-"), sep=":")
            ) %>%
            remove_rownames() %>%
            column_to_rownames("feature")

print(
    paste(
        "Number of differentially bound sites:", nrow(db_sites)
    )
)
print(head(db_sites))

for (i in c("bams", "peaks", "summits", "beds", "read_cnt", "occupancy")){
    metadata[[i]] <- NULL
}

export_plots(db_sites, metadata, args)

print("Exporting differentially bound sites")
export_data(
    data=db_sites,
    location=paste(args$output, "diff_rgns.tsv", sep="_"),
)

row_metadata <- db_sites %>%
                select(
                    "log2FoldChange", "padj",
                    paste0(args$condition1, ".mean"),
                    paste0(args$condition2, ".mean")
                ) %>%
                filter(.$padj <= args$padj) %>%
                arrange(desc(log2FoldChange))
col_metadata <- metadata %>%
                mutate_at(colnames(.), as.vector)
print(head(row_metadata))
print(head(col_metadata))

print(
    paste(
        "Filtering normalized read counts matrix to include",
        "only differential binding sites with padj <= ", args$padj
    )
)
norm_counts_mat <- as.matrix(db_sites[as.vector(rownames(row_metadata)), rownames(col_metadata)])
print("Size of the normalized read counts matrix after filtering")
print(dim(norm_counts_mat))
print(head(norm_counts_mat))

print("Exporting normalized read counts to GCT format")
export_gct(
    counts_mat=norm_counts_mat,
    row_metadata=row_metadata,                                        # includes features as row names
    col_metadata=col_metadata,                                        # includes samples as row names
    location=paste(args$output, "read_cnts.gct", sep="_")
)