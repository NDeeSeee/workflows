#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(modules))
suppressMessages(library(argparse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))


export_all_clustering_plots <- function(seurat_data, suffix, args){
    Idents(seurat_data) <- "new.ident"                                                               # safety measure
    downsampled_to <- analyses$get_min_ident_size(SplitObject(seurat_data, split.by="new.ident"))    # need to split it for consistency
    downsampled_data <- subset(seurat_data, downsample=downsampled_to)
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title=paste("Clustered UMAP projected PCA of GEX datasets. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste("gex_res", current_resolution, sep="."),
            label=TRUE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "gex_umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        graphics$silhouette_plot(
            data=seurat_data,
            reduction="pca",
            dims=args$gexndim,
            downsample=500,
            plot_title=paste("Silhouette scores per cell of downsampled GEX datasets. Max 500 cells per cluster. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste("gex_res", current_resolution, sep="."),
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "gex_silh_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
            graphics$dim_plot(
                data=seurat_data,
                reduction="rnaumap",
                plot_title=paste("Split by identity clustered UMAP projected PCA of GEX datasets. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("gex_res", current_resolution, sep="."),
                split_by="new.ident",
                label=TRUE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "gex_umap_spl_by_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cluster split by identity composition plot of GEX datasets.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("gex_res", current_resolution, sep="."),
                split_by="new.ident",
                x_label="Identity",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "gex_comp_gr_by_clst_spl_by_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by identity split by cluster composition plot of GEX datasets.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution),
                legend_title="Identity",
                group_by="new.ident",
                split_by=paste("gex_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "gex_comp_gr_by_idnt_spl_by_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
        if (seurat_data@meta.data$new.ident != seurat_data@meta.data$condition){
            graphics$dim_plot(
                data=seurat_data,
                reduction="rnaumap",
                plot_title=paste("Split by grouping condition clustered UMAP projected PCA of GEX datasets. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("gex_res", current_resolution, sep="."),
                split_by="condition",
                label=TRUE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "gex_umap_spl_by_cond_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cluster split by condition composition plot of GEX datasets.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("gex_res", current_resolution, sep="."),
                split_by="condition",
                x_label="Condition",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "gex_comp_gr_by_clst_spl_by_cond_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by condition split by cluster composition plot of GEX datasets.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution
                ),
                legend_title="Condition",
                group_by="condition",
                split_by=paste("gex_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "gex_comp_gr_by_cond_spl_by_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
        if ("Phase" %in% colnames(seurat_data@meta.data)){
            graphics$dim_plot(
                data=seurat_data,
                reduction="rnaumap",
                plot_title=paste("Split by cell cycle phase clustered UMAP projected PCA of GEX datasets. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("gex_res", current_resolution, sep="."),
                split_by="Phase",
                label=TRUE,
                label_color="black",
                alpha=0.5,
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "gex_umap_spl_by_ph_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cell cycle phase split by identity composition plot of GEX datasets.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution),
                legend_title="Phase",
                group_by="Phase",
                split_by="new.ident",
                x_label="Identity",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "gex_comp_gr_by_ph_spl_by_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cell cycle phase split by cluster composition plot of GEX datasets.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution),
                legend_title="Phase",
                group_by="Phase",
                split_by=paste("gex_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "gex_comp_gr_by_ph_spl_by_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
    }
    rm(downsampled_data)
    gc(verbose=FALSE)
}


export_all_expression_plots <- function(seurat_data, suffix, args) {
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        Idents(seurat_data) <- paste("gex_res", current_resolution, sep=".")
        graphics$dot_plot(
            data=seurat_data,
            features=args$genes,
            plot_title=paste("Log normalized scaled average gene expression per cluster. Resolution", current_resolution),
            x_label="Genes",
            y_label="Clusters",
            cluster_idents=FALSE,
            rootname=paste(args$output, suffix, "expr_avg_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (length(args$genes) > 0){
            for (i in 1:length(args$genes)) {
                current_gene <- args$genes[i]
                graphics$feature_plot(
                    data=seurat_data,
                    features=current_gene,
                    labels=current_gene,
                    reduction="rnaumap",
                    plot_title=paste("Log normalized gene expression per cell. Resolution", current_resolution),
                    label=TRUE,
                    order=TRUE,
                    max_cutoff="q99",  # to prevent cells with overexpressed gene from distoring the color bar
                    combine_guides="keep",
                    width=800,
                    height=800,
                    rootname=paste(args$output, suffix, "expr_per_cell_res", current_resolution, current_gene, sep="_"),
                    pdf=args$pdf
                )
                graphics$vln_plot(
                    data=seurat_data,
                    features=current_gene,
                    labels=current_gene,
                    plot_title=paste("Log normalized gene expression density per cluster. Resolution", current_resolution),
                    legend_title="Cluster",
                    log=TRUE,
                    pt_size=0,
                    combine_guides="collect",
                    width=800,
                    height=600,
                    palette_colors=graphics$D40_COLORS,
                    rootname=paste(args$output, suffix, "expr_dnst_res", current_resolution, current_gene, sep="_"),
                    pdf=args$pdf
                )
            }
        }
    }
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
}


get_args <- function(){
    parser <- ArgumentParser(description="Seurat GEX Cluster Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file",
            "can be produced by sc_gex_reduce.R script and must include",
            "dimensional reduction information stored in the 'pca' and",
            "'rnaumap' slots. GEX information should be stored in the RNA",
            "assay."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--gexndim",
        help=paste(
            "Dimensionality used when constructing nearest-neighbor graph before",
            "clustering (from 1 to 50). If single value N is provided, use from 1 to N",
            "dimensions. If multiple values are provided, subset to only selected",
            "dimensions. The value will be automatically adjusted based on the",
            "dimensional reduction information stored in the 'pca' slot of the loaded",
            "Seurat object",
            "Default: from 1 to 10"
        ),
        type="integer", default=10, nargs="*"
    )
    parser$add_argument(
        "--ametric",
        help=paste(
            "Distance metric used when constructing nearest-neighbor graph before",
            "clustering.",
            "Default: euclidean"
        ),
        type="character", default="euclidean",
        choices=c(
            "euclidean", "cosine", "manhattan", "hamming"
        )
    )
    parser$add_argument(
        "--resolution",
        help=paste(
            "Clustering resolution applied to the constructed nearest-neighbor graph.",
            "Can be set as an array.",
            "Default: 0.3, 0.5, 1.0"
        ),
        type="double", default=c(0.3, 0.5, 1.0), nargs="*"
    )
    parser$add_argument(
        "--genes",
        help=paste(
            "Genes of interest to evaluate expression.",
            "Default: None"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--gexpttv",
        help=paste(
            "Identify putative gene markers for all clusters and all resolutions.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--gexlogfc",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "on average have log fold change difference in expression between every",
            "tested pair of clusters not lower than this value. Ignored if --gexpttv",
            "is not set.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--gexminpct",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "are detected in not lower than this fraction of cells in either of the",
            "two tested clusters. Ignored if --gexpttv is not set.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--gexonlypos",
        help=paste(
            "For putative gene markers identification return only positive markers.",
            "Ignored if --gexpttv is not set.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--gextestuse",
        help=paste(
            "Statistical test to use for putative gene markers identification.",
            "Ignored if --gexpttv is not set.",
            "Default: wilcox"
        ),
        type="character", default="wilcox",
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
    )
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--verbose",
        help="Print debug information. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--h5seurat",
        help="Save Seurat data to h5seurat file. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--cbbuild",
        help="Export results to UCSC Cell Browser. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./seurat",
        type="character", default="./seurat"
    )
    parser$add_argument(
        "--cpus",
        help="Number of cores/cpus to use. Default: 1",
        type="integer", default=1
    )
    parser$add_argument(
        "--memory",
        help=paste(
            "Maximum memory in GB allowed to be shared between the workers",
            "when using multiple --cpus.",
            "Default: 32"
        ),
        type="integer", default=32
    )
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}

args <- get_args()

print("Input parameters")
print(args)
if (length(args$gexndim) == 1) {
    print("Adjusting --gexndim parameter as only a single value was provided")
    args$gexndim <- c(1:args$gexndim[1])
    print(paste("--gexndim was adjusted to", paste(args$gexndim, collapse=", ")))
}

print(
    paste(
        "Setting parallelization to", args$cpus, "cores, and", args$memory,
        "GB of memory allowed to be shared between the processes"
    )
)
prod$parallel(args)

print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)
print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"
debug$print_info(seurat_data, args)

if (!all(c("pca", "rnaumap") %in% names(seurat_data@reductions))){
    print("Loaded Seurat object doesn't have 'pca' and/or 'rnaumap' reduction(s). Exiting.")
    quit(save="no", status=1, runLast=FALSE)
}

if (!is.null(args$genes)){
    print("Adjusting genes of interest to include only those that are present in the loaded Seurat object")
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]     # with RNA assay set as default the rownames should be genes
    print(args$genes)
}

print(paste("Clustering GEX data using", paste(args$gexndim, collapse=", "), "principal components"))
seurat_data <- analyses$gex_cluster(
    seurat_data=seurat_data,
    graph_name="gex",                          # will be used on all the plot generating functions
    args=args
)
debug$print_info(seurat_data, args)

export_all_clustering_plots(
    seurat_data=seurat_data,
    suffix="clst",
    args=args
)

if(args$cbbuild){
    print("Exporting UCSC Cellbrowser data")
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        features=args$genes,                                   # can be NULL
        resolution=args$resolution,
        resolution_prefix="gex_res",
        rootname=paste(args$output, "_cellbrowser", sep=""),
    )
}

if (!is.null(args$genes) || args$gexpttv) {
    print("Normalizing counts in RNA assay before evaluating genes expression or identifying putative gene markers")
    DefaultAssay(seurat_data) <- "RNA"
    seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
    if (!is.null(args$genes)){
        print("Generating plot to evaluate genes expression")
        export_all_expression_plots(seurat_data=seurat_data, suffix="clst", args=args)
    }
    if(args$gexpttv){
        print("Identifying putative gene markers for all clusters and all resolutions")
        all_putative_markers <- analyses$gex_putative_gene_markers(
            seurat_data=seurat_data,
            resolution_prefix="gex_res",
            args=args
        )
        io$export_data(
            all_putative_markers,
            paste(args$output, "_clst_pttv_gene_markers.tsv", sep="")
        )
    }
}

DefaultAssay(seurat_data) <- "RNA"                                                         # better to stick to RNA assay by default https://www.biostars.org/p/395951/#395954 
io$export_rds(seurat_data, paste(args$output, "_clst_data.rds", sep=""))
if(args$h5seurat){
    io$export_h5seurat(seurat_data, paste(args$output, "_clst_data.h5seurat", sep=""))
}