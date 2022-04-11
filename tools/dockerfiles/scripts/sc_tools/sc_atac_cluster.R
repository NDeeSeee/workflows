#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(modules))
suppressMessages(library(argparse))
suppressMessages(library(GenomicRanges))

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
            reduction="atacumap",
            plot_title=paste("Clustered UMAP projected LSI of ATAC datasets. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste("atac_res", current_resolution, sep="."),
            label=TRUE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "atac_umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        graphics$silhouette_plot(
            data=seurat_data,
            reduction="atac_lsi",
            dims=args$atacndim,
            downsample=500,
            plot_title=paste("Silhouette scores per cell of downsampled ATAC datasets. Max 500 cells per cluster. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste("atac_res", current_resolution, sep="."),
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "atac_silh_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
            graphics$dim_plot(
                data=seurat_data,
                reduction="atacumap",
                plot_title=paste("Split by identity clustered UMAP projected LSI of ATAC datasets. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("atac_res", current_resolution, sep="."),
                split_by="new.ident",
                label=TRUE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "atac_umap_spl_by_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cluster split by identity composition plot of ATAC datasets.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("atac_res", current_resolution, sep="."),
                split_by="new.ident",
                x_label="Identity",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "atac_comp_gr_by_clst_spl_by_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by identity split by cluster composition plot of ATAC datasets.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution),
                legend_title="Identity",
                group_by="new.ident",
                split_by=paste("atac_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "atac_comp_gr_by_idnt_spl_by_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
        if (seurat_data@meta.data$new.ident != seurat_data@meta.data$condition){
            graphics$dim_plot(
                data=seurat_data,
                reduction="atacumap",
                plot_title=paste("Split by grouping condition clustered UMAP projected LSI of ATAC datasets. Resolution", current_resolution),
                legend_title="Cluster",
                group_by=paste("atac_res", current_resolution, sep="."),
                split_by="condition",
                label=TRUE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "atac_umap_spl_by_cond_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by cluster split by condition composition plot of ATAC datasets.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("atac_res", current_resolution, sep="."),
                split_by="condition",
                x_label="Condition",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "atac_comp_gr_by_clst_spl_by_cond_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Grouped by condition split by cluster composition plot of ATAC datasets.",
                    "Downsampled to", downsampled_to, "cells per dataset.",
                    "Resolution", current_resolution
                ),
                legend_title="Condition",
                group_by="condition",
                split_by=paste("atac_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "atac_comp_gr_by_cond_spl_by_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
    }
    rm(downsampled_data)
    gc(verbose=FALSE)
}


export_all_coverage_plots <- function(seurat_data, suffix, show_expression, args) {
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"                                          # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                           # safety measure

    genome_annotation <- Annotation(seurat_data)                                               # safety measure to build the coverage plot
    if( !("gene_biotype" %in% base::colnames(GenomicRanges::mcols(genome_annotation))) ){
        print("Genome annotation doesn't have 'gene_biotype' column. Adding NA")
        genome_annotation$gene_biotype <- NA
        Annotation(seurat_data) <- genome_annotation
    }

    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        for (i in 1:length(args$genes)) {
            current_gene <- args$genes[i]
            graphics$coverage_plot(
                data=seurat_data,
                assay="ATAC",
                region=current_gene,
                group_by=paste("atac_res", current_resolution, sep="."),
                plot_title=paste(
                    "Tn5 insertion frequency around", current_gene, "gene grouped by cluster for all datasets.",
                    "Resolution", current_resolution
                ),
                idents=NULL,                                                   # to include all values from the default "new.ident" column
                cells=colnames(seurat_data),                                   # limit to only those cells that are in out seurat_data
                features=if(show_expression) current_gene else NULL,           # otherwise will fail if features are set but RNA assay is not present
                expression_assay="RNA",
                expression_slot="data",                                        # use scaled counts
                extend_upstream=2500,
                extend_downstream=2500,
                show_annotation=TRUE,
                show_peaks=TRUE,
                palette_colors=graphics$D40_COLORS,
                rootname=paste(args$output, suffix, "cvrg_res", current_resolution, current_gene, sep="_"),
                pdf=args$pdf
            )
        }
    }
}


get_args <- function(){
    parser <- ArgumentParser(description="Seurat ATAC Cluster Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file",
            "can be produced by sc_atac_reduce.R script and must include",
            "dimensional reduction information stored in the 'atac_lsi' and",
            "'atacumap' slots. It is mandatory to have chromatin accessibility",
            "information stored in the ATAC assay."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--atacndim",
        help=paste(
            "Dimensionality used when constructing nearest-neighbor graph before",
            "clustering (from 1 to 50). If single value N is provided, use from 2 to N",
            "dimensions. If multiple values are provided, subset to only selected",
            "dimensions.",
            "Default: from 2 to 10"
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
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment used in the",
            "construction of Seurat object loaded with --query parameter.",
            "File should be saved in TSV format with tbi-index file."
        ),
        type="character"
    )
    parser$add_argument(
        "--genes",
        help=paste(
            "Genes of interest to build Tn5 insertion frequency plots. If loaded Seurat object",
            "includes RNA assay it will be additionally shown on the side of the plots.",
            "Ignored if --fragments is not provided.",
            "Default: None"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--diffpeaks",
        help=paste(
            "Identify differentially accessible peaks for all clusters and all resolutions.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--ataclogfc",
        help=paste(
            "For differentially accessible peaks identification include only those peaks that",
            "on average have log fold change difference in the chromatin accessibility between",
            "every tested pair of clusters not lower than this value. Ignored if --diffpeaks",
            "is not set.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--atacminpct",
        help=paste(
            "For differentially accessible peaks identification include only those peaks that",
            "are detected in not lower than this fraction of cells in either of the two tested",
            "clusters. Ignored if --diffpeaks is not set.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--atactestuse",
        help=paste(
            "Statistical test to use for differentially accessible peaks identification.",
            "Ignored if --diffpeaks is not set.",
            "Default: LR"
        ),
        type="character", default="LR",
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
if (length(args$atacndim) == 1) {
    print("Adjusting --atacndim parameter as only a single value was provided")
    args$atacndim <- c(2:args$atacndim[1])                                              # skipping the first LSI component
    print(paste("--atacndim was adjusted to", paste(args$atacndim, collapse=", ")))
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
print("Setting default assay to ATAC")
DefaultAssay(seurat_data) <- "ATAC"
debug$print_info(seurat_data, args)

if (!all(c("atac_lsi", "atacumap") %in% names(seurat_data@reductions))){
    print("Loaded Seurat object doesn't have 'atac_lsi' and/or 'atacumap' reduction(s). Exiting.")
    quit(save="no", status=1, runLast=FALSE)
}

if (!is.null(args$fragments)){
    print(paste("Loading fragments data from", args$fragments))
    seurat_data <- io$replace_fragments(args$fragments, seurat_data)
    debug$print_info(seurat_data, args)
}

print(paste("Clustering ATAC data using", paste(args$atacndim, collapse=", "), "dimensions"))
seurat_data <- analyses$add_clusters(
    seurat_data=seurat_data,
    assay="ATAC",
    graph_name="atac",                          # will be used on all the plot generating functions
    reduction="atac_lsi",
    dims=args$atacndim,
    cluster_algorithm=3,                        # SLM algorithm
    args=args
)
debug$print_info(seurat_data, args)

export_all_clustering_plots(
    seurat_data=seurat_data,
    suffix="clst",
    args=args
)

regions <- NULL
if (!is.null(args$genes)){
    print("Adjusting genes of interest to include only those that are present in the loaded Seurat object")
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(Annotation(seurat_data)$gene_name))]
    regions <- sapply(
        args$genes,
        function(gene)
        GRangesToString(
            grange=LookupGeneCoords(
                seurat_data,
                gene=gene,
                assay="ATAC"
            ),
            sep=c("-", "-")
        )
    )
    print(regions)
}

if(args$cbbuild){
    print("Exporting UCSC Cellbrowser data")
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="ATAC",
        slot="counts",
        features=regions,                                      # use regions that correspond to our genes of interest
        rootname=paste(args$output, "_cellbrowser", sep=""),
    )
}

if (!is.null(args$genes) && !is.null(args$fragments)){
    show_expression <- FALSE
    if ("RNA" %in% names(seurat_data@assays)){
        print("Normalizing counts in RNA assay to show average gene expression alongside the coverage plots")
        DefaultAssay(seurat_data) <- "RNA"                         # for now we will always use RNA even if SCT can be present
        seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
        DefaultAssay(seurat_data) <- "ATAC"                        # safety measure
        show_expression <- TRUE
    }
    export_all_coverage_plots(
        seurat_data=seurat_data,
        suffix="clst",
        show_expression=show_expression,
        args=args
    )
}

if (args$diffpeaks){
    print("Running differential accesibility analysis")
}

DefaultAssay(seurat_data) <- "ATAC"
io$export_rds(seurat_data, paste(args$output, "_clst_data.rds", sep=""))
if(args$h5seurat){
    io$export_h5seurat(seurat_data, paste(args$output, "_clst_data.h5seurat", sep=""))
}