#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(modules))
suppressMessages(library(forcats))
suppressMessages(library(argparse))
suppressMessages(library(tidyselect))
suppressMessages(library(GenomicRanges))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))


export_all_clustering_plots <- function(seurat_data, args){
    Idents(seurat_data) <- "new.ident"                                                               # safety measure
    downsampled_to <- analyses$get_min_ident_size(SplitObject(seurat_data, split.by="new.ident"))    # need to split it for consistency
    print(paste("Downsampling datasets to", downsampled_to, "cells per sample"))
    downsampled_data <- subset(seurat_data, downsample=downsampled_to)
    if ("Phase" %in% colnames(seurat_data@meta.data)){
        if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Composition plot,",
                    "colored by cell cycle phase,",
                    "split by dataset, downsampled"
                ),
                legend_title="Phase",
                group_by="Phase",
                split_by="new.ident",
                x_label="Dataset",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_ph_spl_idnt", sep="_"),
                pdf=args$pdf
            )
        }
    }
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        graphics$dim_plot(
            data=seurat_data,
            reduction="wnnumap",
            plot_title=paste(
                "UMAP, colored by cluster,",
                "resolution", current_resolution
            ),
            legend_title="Cluster",
            group_by=paste("wsnn_res", current_resolution, sep="."),
            label=TRUE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
            graphics$dim_plot(
                data=seurat_data,
                reduction="wnnumap",
                plot_title=paste(
                    "UMAP, colored by cluster,",
                    "split by dataset,",
                    "resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("wsnn_res", current_resolution, sep="."),
                split_by="new.ident",
                label=TRUE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Composition plot, colored by cluster,",
                    "split by dataset, downsampled,",
                    "resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("wsnn_res", current_resolution, sep="."),
                split_by="new.ident",
                x_label="Dataset",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_clst_spl_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Composition plot, colored by dataset,",
                    "split by cluster, downsampled,",
                    "resolution", current_resolution
                ),
                legend_title="Dataset",
                group_by="new.ident",
                split_by=paste("wsnn_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_idnt_spl_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
        if (
            all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
            length(unique(as.vector(as.character(seurat_data@meta.data$condition)))) > 1
        ){
            graphics$dim_plot(
                data=seurat_data,
                reduction="wnnumap",
                plot_title=paste(
                    "UMAP, colored by cluster,",
                    "split by grouping condition,",
                    "resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("wsnn_res", current_resolution, sep="."),
                split_by="condition",
                label=TRUE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_cnd_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Composition plot, colored by cluster,",
                    "split by grouping condition, downsampled,",
                    "resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("wsnn_res", current_resolution, sep="."),
                split_by="condition",
                x_label="Condition",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_clst_spl_cnd_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Composition plot,",
                    "colored by grouping condition,",
                    "split by cluster, downsampled,",
                    "resolution", current_resolution
                ),
                legend_title="Condition",
                group_by="condition",
                split_by=paste("wsnn_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_cnd_spl_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
        if ("Phase" %in% colnames(seurat_data@meta.data)){
            graphics$dim_plot(
                data=seurat_data,
                reduction="wnnumap",
                plot_title=paste(
                    "UMAP, colored by cluster,",
                    "split by cell cycle phase,",
                    "resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("wsnn_res", current_resolution, sep="."),
                split_by="Phase",
                label=TRUE,
                label_color="black",
                alpha=0.5,
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_spl_ph_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_data,
                plot_title=paste(
                    "Composition plot,",
                    "colored by cell cycle phase,",
                    "split by cluster, downsampled,",
                    "resolution", current_resolution
                ),
                legend_title="Phase",
                group_by="Phase",
                split_by=paste("wsnn_res", current_resolution, sep="."),
                x_label="Cluster",
                y_label="Cells percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_ph_spl_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
        }
    }
    rm(downsampled_data)
    gc(verbose=FALSE)
}


export_all_expression_plots <- function(seurat_data, args) {
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    if (!is.null(args$genes) && length(args$genes) > 0){
        for (i in 1:length(args$genes)){
            current_gene <- args$genes[i]
            graphics$feature_plot(
                data=seurat_data,
                features=current_gene,
                labels=current_gene,
                reduction="wnnumap",
                plot_title="UMAP, gene expression",
                label=FALSE,
                order=TRUE,
                max_cutoff="q99",  # to prevent cells with overexpressed gene from distoring the color bar
                combine_guides="keep",
                width=800,
                height=800,
                theme=args$theme,
                rootname=paste(args$output, "xpr_per_cell", current_gene, sep="_"),
                pdf=args$pdf
            )
            graphics$expression_density_plot(
                data=seurat_data,
                features=current_gene,
                reduction="wnnumap",
                plot_title="UMAP, gene expression density",
                joint=FALSE,
                width=800,
                height=800,
                theme=args$theme,
                rootname=paste(args$output, "xpr_per_cell_sgnl", current_gene, sep="_"),
                pdf=args$pdf
            )
        }
    }
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        Idents(seurat_data) <- paste("wsnn_res", current_resolution, sep=".")
        graphics$dot_plot(
            data=seurat_data,
            features=args$genes,
            plot_title=paste(
                "Gene expression dot plot,",
                "resolution", current_resolution
            ),
            x_label="Genes",
            y_label="Clusters",
            cluster_idents=FALSE,
            theme=args$theme,
            rootname=paste(args$output, "xpr_avg_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (!is.null(args$genes) && length(args$genes) > 0){
            for (i in 1:length(args$genes)){
                current_gene <- args$genes[i]
                graphics$vln_plot(
                    data=seurat_data,
                    features=current_gene,
                    labels=current_gene,
                    plot_title=paste(
                        "Gene expression violin plot,",
                        "resolution", current_resolution
                    ),
                    legend_title="Cluster",
                    log=TRUE,
                    pt_size=0,
                    combine_guides="collect",
                    width=800,
                    height=600,
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "xpr_dnst_res", current_resolution, current_gene, sep="_"),
                    pdf=args$pdf
                )
            }
        }
    }
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
}


export_all_coverage_plots <- function(seurat_data, args) {
    SeuratObject::DefaultAssay(seurat_data) <- "ATAC"                                          # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                                           # safety measure

    genome_annotation <- Annotation(seurat_data)                                               # safety measure to build the coverage plot
    if( !("gene_biotype" %in% base::colnames(GenomicRanges::mcols(genome_annotation))) ){
        print("Genome annotation doesn't have 'gene_biotype' column. Adding NA")
        genome_annotation$gene_biotype <- NA
    }
    if( !("tx_id" %in% base::colnames(GenomicRanges::mcols(genome_annotation))) ){             # https://github.com/stuart-lab/signac/issues/1159
        print("Genome annotation doesn't have 'tx_id' column. Adding from 'transcript_id'")
        genome_annotation$tx_id <- genome_annotation$transcript_id
    }
    Annotation(seurat_data) <- genome_annotation

    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        if (!is.null(args$genes) && length(args$genes) > 0){
            for (i in 1:length(args$genes)){
                current_gene <- args$genes[i]
                graphics$coverage_plot(
                    data=seurat_data,
                    assay="ATAC",
                    region=current_gene,
                    group_by=paste("wsnn_res", current_resolution, sep="."),
                    plot_title=paste(
                        "Fragments coverage,",
                        current_gene, "gene,",
                        "resolution", current_resolution
                    ),
                    idents=NULL,                                                   # to include all values from the default "new.ident" column
                    cells=colnames(seurat_data),                                   # limit to only those cells that are in out seurat_data
                    features=current_gene,
                    expression_assay="RNA",
                    expression_slot="data",                                        # use scaled counts
                    extend_upstream=2500,
                    extend_downstream=2500,
                    show_annotation=TRUE,
                    show_peaks=TRUE,
                    show_tile=TRUE,
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "cvrg_res", current_resolution, current_gene, sep="_"),
                    pdf=args$pdf
                )
            }
        }
    }
}


export_heatmaps <- function(seurat_data, markers, args){
    DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        cluster_column <- paste("wsnn_res", current_resolution, sep=".")
        clusters_order <- levels(seurat_data@meta.data[[cluster_column]])
        if (is.null(clusters_order)){
            clusters_order <- unique(seurat_data@meta.data[[cluster_column]])
        }

        filtered_markers <- markers %>%
                            dplyr::filter(.$resolution==current_resolution) %>%                 # shouldn't fail even if resolution is not present
                            dplyr::select(-c("resolution")) %>%
                            dplyr::mutate(
                                cluster=base::factor(
                                    cluster,
                                    levels=clusters_order
                                )
                            ) %>%
                            dplyr::arrange(cluster) %>%                                         # to have upregulated genes on the diagonal
                            dplyr::filter(.$p_val_adj <= 0.05) %>%                              # to have only significant gene markers
                            dplyr::filter(.$pct.1 >= 0.1) %>%                                   # to have at least 10% of cells expressing this gene
                            dplyr::group_by(feature) %>%
                            dplyr::arrange(desc(pct.1), .by_group=TRUE) %>%                     # sort all duplicated features by desc(pct.1)
                            dplyr::slice_head(n=1) %>%                                          # choose the feature with the highest pct.1
                            dplyr::ungroup() %>%
                            dplyr::group_by(cluster) %>%
                            dplyr::arrange(p_val_adj, desc(avg_log2FC), .by_group=TRUE) %>%
                            dplyr::group_modify(~ .x %>%
                                slice_head(n=analyses$get_fraction(.x, 0.25))                   # take 25% of the features
                            ) %>%
                            dplyr::arrange(desc(avg_log2FC), .by_group=TRUE) %>%
                            dplyr::ungroup()

        if (nrow(filtered_markers) > 0){                                                            # in case we don't have any markers for specific resolution
            column_annotations <- c(cluster_column)
            if (
                all(as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))) &&
                length(unique(as.vector(as.character(seurat_data@meta.data$condition)))) > 1
            ){
                column_annotations <- c(column_annotations, "condition")                           # several conditions found
            }
            if (length(unique(as.vector(as.character(seurat_data@meta.data$new.ident)))) > 1){
                column_annotations <- c(column_annotations, "new.ident")                           # several datasets found
            }
            graphics$feature_heatmap(                                                              # install.packages("magick") for better rasterization
                data=seurat_data,
                assay="RNA",
                slot="data",
                features=filtered_markers$feature,
                highlight_features=if(!is.null(args$genes) && length(args$genes) > 0) args$genes else NULL,
                show_rownames=FALSE,
                scale_to_max=FALSE,
                scale="row",                                                                       # will calculate z-score
                heatmap_colors=c("darkblue", "black", "yellow"),
                group_by=unique(column_annotations),                                               # safety measure for possible duplicates
                palette_colors=graphics$D40_COLORS,
                plot_title=paste(
                    "Gene expression heatmap, resolution", current_resolution
                ),
                rootname=paste(args$output, "xpr_htmp_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            io$export_data(
                filtered_markers,
                paste(args$output, "_xpr_htmp_res_", current_resolution, ".tsv", sep="")
            )
        }
    }
    Idents(seurat_data) <- "new.ident"                            # safety measure
}


get_args <- function(){
    parser <- ArgumentParser(description="Single-cell WNN Cluster Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include",
            "genes expression and chromatin accessibility information stored in the RNA",
            "and ATAC assays correspondingly. Additionally, 'pca', 'rnaumap', 'atac_lsi'",
            "and 'atacumap' dimensionality reductions should be present."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--rnadimensions",
        help=paste(
            "Dimensionality from the 'pca' reduction to use when constructing weighted",
            "nearest-neighbor graph before clustering (from 1 to 50).",
            "Default: 10"
        ),
        type="integer", default=10
    )
    parser$add_argument(
        "--atacdimensions",
        help=paste(
            "Dimensionality from the 'atac_lsi' reduction to use when constructing weighted",
            "nearest-neighbor graph before clustering (from 2 to 50). First LSI component is",
            "always excluded unless the provided RDS file consists of multiple datasets",
            "where ATAC assay were integrated with Harmony.",
            "Default: 10"
        ),
        type="integer", default=10
    )
    parser$add_argument(
        "--algorithm",
        help=paste(
            "Algorithm for modularity optimization when running clustering.",
            "Default: louvain"
        ),
        type="character", default="slm",
        choices=c(
            "louvain", "mult-louvain", "slm", "leiden"
        )
    )
    parser$add_argument(
        "--uspread",
        help=paste(
            "The effective scale of embedded points on UMAP. In combination with '--mindist'",
            "it determines how clustered/clumped the embedded points are.",
            "Default: 1"
        ),
        type="double", default=1
    )
    parser$add_argument(
        "--umindist",
        help=paste(
            "Controls how tightly the embedding is allowed compress points together on UMAP.",
            "Larger values ensure embedded points are moreevenly distributed, while smaller",
            "values allow the algorithm to optimise more accurately with regard to local structure.",
            "Sensible values are in the range 0.001 to 0.5.",
            "Default:  0.3"
        ),
        type="double", default=0.3
    )
    parser$add_argument(
        "--uneighbors",
        help=paste(
            "Determines the number of neighboring points used in UMAP. Larger values will result",
            "in more global structure being preserved at the loss of detailed local structure.",
            "In general this parameter should often be in the range 5 to 50.",
            "Default: 30"
        ),
        type="integer", default=30
    )
    parser$add_argument(
        "--umetric",
        help=paste(
            "The metric to use to compute distances in high dimensional space for UMAP.",
            "Default: cosine"
        ),
        type="character", default="cosine",
        choices=c(
            "euclidean", "manhattan", "chebyshev", "minkowski", "canberra", "braycurtis",
            "mahalanobis", "wminkowski", "seuclidean", "cosine", "correlation", "haversine",
            "hamming", "jaccard", "dice", "russelrao", "kulsinski", "ll_dirichlet", "hellinger",
            "rogerstanimoto", "sokalmichener", "sokalsneath", "yule"
        )
    )
    # The default method for RunUMAP has changed from calling Python UMAP via reticulate to
    # the R-native UWOT using the cosine metric. To use Python UMAP via reticulate, set
    # umap.method to 'umap-learn' and metric to 'correlation'
    parser$add_argument(
        "--umethod",
        help=paste(
            "UMAP implementation to run. If set to 'umap-learn' use --umetric 'correlation'",
            "Default: uwot"
        ),
        type="character", default="uwot",
        choices=c("uwot", "uwot-learn", "umap-learn")
    )
    parser$add_argument(
        "--resolution",
        help=paste(
            "Clustering resolution applied to the constructed weighted nearest-neighbor",
            "graph. Can be set as an array but only the first item from the list will",
            "be used for cluster labels and gene/peak markers in the UCSC Cell Browser",
            "when running with --cbbuild and --diffgenes/--diffpeaks parameters.",
            "Default: 0.3, 0.5, 1.0"
        ),
        type="double", default=c(0.3, 0.5, 1.0), nargs="*"
    )
    parser$add_argument(
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment used in the loaded Seurat",
            "object. File should be saved in TSV format with tbi-index file."
        ),
        type="character"
    )
    parser$add_argument(
        "--genes",
        help=paste(
            "Genes of interest to build gene expression and Tn5 insertion frequency plots",
            "for the nearest peaks. If '--fragments' is not provided only gene expression",
            "plots will be built.",
            "Default: None"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--diffgenes",
        help=paste(
            "Identify differentially expressed genes (putative gene markers) between each",
            "pair of clusters for all resolutions.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--diffpeaks",
        help=paste(
            "Identify differentially accessible peaks between each pair of clusters for all resolutions.",
            "Default: false"
        ),
        action="store_true"
    )

    parser$add_argument(
        "--rnalogfc",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "on average have log fold change difference in expression between every",
            "tested pair of clusters not lower than this value. Ignored if '--diffgenes'",
            "is not set.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--rnaminpct",
        help=paste(
            "For putative gene markers identification include only those genes that",
            "are detected in not lower than this fraction of cells in either of the",
            "two tested clusters. Ignored if '--diffgenes' is not set.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--rnaonlypos",
        help=paste(
            "For putative gene markers identification return only positive markers.",
            "Ignored if '--diffgenes' is not set.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--rnatestuse",
        help=paste(
            "Statistical test to use for putative gene markers identification.",
            "Ignored if '--diffgenes' is not set.",
            "Default: wilcox"
        ),
        type="character", default="wilcox",
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
    )
    parser$add_argument(
        "--ataclogfc",
        help=paste(
            "For differentially accessible peaks identification include only those peaks that",
            "on average have log fold change difference in the chromatin accessibility between",
            "every tested pair of clusters not lower than this value. Ignored if '--diffpeaks'",
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
            "clusters. Ignored if '--diffpeaks' is not set.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--atactestuse",
        help=paste(
            "Statistical test to use for differentially accessible peaks identification.",
            "Ignored if '--diffpeaks' is not set.",
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
        "--h5ad",
        help="Save Seurat data to h5ad file. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--cbbuild",
        help="Export results to UCSC Cell Browser. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--scope",
        help=paste(
            "Save Seurat data to SCope compatible loom file. Only",
            "not normalized raw counts from the RNA assay will be",
            "saved. Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./sc",
        type="character", default="./sc"
    )
    parser$add_argument(
        "--theme",
        help=paste(
            "Color theme for all generated plots.",
            "Default: classic"
        ),
        type="character", default="classic",
        choices=c("gray", "bw", "linedraw", "light", "dark", "minimal", "classic", "void")
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
DefaultAssay(seurat_data) <- "RNA"                                                               # for consistency with other scripts
debug$print_info(seurat_data, args)

if (!all(c("pca", "atac_lsi") %in% names(seurat_data@reductions))){
    print("Loaded Seurat object doesn't have 'pca' and/or 'atac_lsi' reduction(s). Exiting.")
    quit(save="no", status=1, runLast=FALSE)
}

print("Adjusting --rnadimensions parameter")
args$rnadimensions <- c(1:args$rnadimensions)
print(paste("--rnadimensions was adjusted to", paste(args$rnadimensions, collapse=", ")))

print("Adjusting --atacdimensions parameter")
if (!is.null(seurat_data@misc$atac_reduce$first_lsi_removed) && seurat_data@misc$atac_reduce$first_lsi_removed){
    args$atacdimensions <- c(1:(max(args$atacdimensions)-1))  # first LSI component has been already removed
} else {
    args$atacdimensions <- c(2:args$atacdimensions)
}
print(paste("--atacdimensions was adjusted to", paste(args$atacdimensions, collapse=", ")))

if (!is.null(args$fragments)){
    print(paste("Loading fragments data from", args$fragments))
    seurat_data <- io$replace_fragments(args$fragments, seurat_data)                             # will change the default assay to ATAC
    debug$print_info(seurat_data, args)
}

print(
    paste(
        "Running weighted nearest-neighbor analysis using", paste(args$rnadimensions, collapse=", "),
        "dimensions from 'pca' and", paste(args$atacdimensions, collapse=", "), "dimensions from 'atac_lsi'",
        "dimensionality reductions."
    )
)
seurat_data <- analyses$add_wnn_clusters(                       # will add 'wnnumap' reduction
    seurat_data=seurat_data,
    graph_name="wsnn",                                          # will be used in all the plot generating functions
    reductions=list("pca", "atac_lsi"),
    dimensions=list(args$rnadimensions, args$atacdimensions),   # should be the same order as reductions
    args=args
)
debug$print_info(seurat_data, args)

export_all_clustering_plots(seurat_data=seurat_data, args=args)

nearest_peaks <- NULL
if (!is.null(args$genes)){
    print("Adjusting genes of interest to include only those that are present in the loaded Seurat object")
    args$genes <- unique(args$genes)
    DefaultAssay(seurat_data) <- "RNA"                                                                       # need it for rownames to return genes
    args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]                 # with RNA assay set as default the rownames should be genes
    DefaultAssay(seurat_data) <- "ATAC"                                                                      # Annotation needs the default assay to be ATAC
    args$genes <- args$genes[args$genes %in% as.vector(as.character(Annotation(seurat_data)$gene_name))]     # just in case check if the same genes are present in Annotation
    all_peaks <- StringToGRanges(rownames(seurat_data), sep=c("-", "-"))                                     # rownames are peaks when default assay is ATAC
    nearest_peaks <- sapply(
        args$genes,
        function(gene)
        GRangesToString(
            all_peaks[
                subjectHits(
                    distanceToNearest(
                        LookupGeneCoords(seurat_data, gene=gene, assay="ATAC"),
                        all_peaks
                    )
                )
            ]
        )
    )
    rm(all_peaks)
    print(nearest_peaks)
}

all_rna_markers <- NULL
if (!is.null(args$genes) || args$diffgenes) {
    print("Normalizing counts in RNA assay")
    DefaultAssay(seurat_data) <- "RNA"
    seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
    if (!is.null(args$genes)){
        print("Generating genes expression plots")
        export_all_expression_plots(seurat_data=seurat_data, args=args)
    }
    if (!is.null(args$genes) && !is.null(args$fragments)){
        print("Generating coverage plots")
        export_all_coverage_plots(seurat_data=seurat_data, args=args)
    }
    if(args$diffgenes){
        print("Identifying differentially expressed genes between each pair of clusters for all resolutions")
        args$logfc <- args$rnalogfc                                  # need the proper names for get_markers_by_res
        args$minpct <- args$rnaminpct
        args$onlypos <- args$rnaonlypos
        args$testuse <- args$rnatestuse
        all_rna_markers <- analyses$get_markers_by_res(              # will change default assay to RNA
            seurat_data=seurat_data,
            assay="RNA",
            resolution_prefix="wsnn_res",
            args=args
        )
        args <- args[names(args) %in% c("logfc", "minpct", "onlypos", "testuse") == FALSE]  # to remove temporary added items
        if (!is.null(all_rna_markers)){
            io$export_data(
                all_rna_markers,
                paste(args$output, "_gene_markers.tsv", sep="")
            )
            export_heatmaps(                                         # will change default assay to RNA
                seurat_data=seurat_data,
                markers=all_rna_markers,
                args=args
            )
        }
    }
}

all_atac_markers <- NULL
if (args$diffpeaks){
    print("Identifying differentially accessible peaks between each pair of clusters for all resolutions")
    DefaultAssay(seurat_data) <- "ATAC"                              # safety measure
    args$logfc <- args$ataclogfc                                     # need the proper names for get_markers_by_res
    args$minpct <- args$atacminpct
    args$onlypos <- FALSE                                            # need to overwrite what was set for RNA
    args$testuse <- args$atactestuse
    all_atac_markers <- analyses$get_markers_by_res(                 # will change default assay to ATAC
        seurat_data=seurat_data,
        assay="ATAC",
        resolution_prefix="wsnn_res",
        latent_vars="nCount_ATAC",                                   # to remove the influence of sequencing depth
        args=args
    )
    args <- args[names(args) %in% c("logfc", "minpct", "onlypos", "testuse") == FALSE]  # to remove temporary added items
    if (!is.null(all_atac_markers)){
        io$export_data(
            all_atac_markers,
            paste(args$output, "_peak_markers.tsv", sep="")
        )
    }
}

if(args$cbbuild){
    print("Exporting RNA and ATAC assays to UCSC Cellbrowser")

    if(!is.null(all_rna_markers)){
        all_rna_markers <- all_rna_markers %>%
                           dplyr::filter(.$resolution==args$resolution[1]) %>%     # won't fail even if resolution is not present
                           dplyr::select(-c("resolution"))
    }
    if(!is.null(all_atac_markers)){
        all_atac_markers <- all_atac_markers %>%
                            dplyr::filter(.$resolution==args$resolution[1]) %>%    # won't fail even if resolution is not present
                            dplyr::select(-c("resolution"))
    }

    print("Reordering reductions to have wnnumap on the first place")              # will be shown first in UCSC Cellbrowser
    reduc_names <- names(seurat_data@reductions)
    ordered_reduc_names <- c("wnnumap", reduc_names[reduc_names!="wnnumap"])       # wnnumap will be added by this time
    seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]
    debug$print_info(seurat_data, args)

    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        short_label="RNA",
        markers=all_rna_markers,                                                   # can be NULL
        label_field=paste0("Clustering (wsnn ", args$resolution[1], ")"),          # always use only the first resolution
        is_nested=TRUE,
        palette_colors=graphics$D40_COLORS,                                        # to have colors correspond to the plots
        rootname=paste(args$output, "_cellbrowser/rna", sep="")
    )

    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="ATAC",
        slot="counts",
        short_label="ATAC",
        markers=all_atac_markers,                                                  # can be NULL
        label_field=paste0("Clustering (wsnn ", args$resolution[1], ")"),          # always use only the first resolution
        is_nested=TRUE,
        palette_colors=graphics$D40_COLORS,                                        # to have colors correspond to the plots
        rootname=paste(args$output, "_cellbrowser/atac", sep="")
    )
}

DefaultAssay(seurat_data) <- "RNA"
print("Exporting results to RDS file")
io$export_rds(seurat_data, paste(args$output, "_data.rds", sep=""))
if(args$h5seurat){
    print("Exporting results to h5seurat file")
    io$export_h5seurat(seurat_data, paste(args$output, "_data.h5seurat", sep=""))
}

if(args$h5ad){
    print("Exporting results to h5ad file")
    io$export_h5ad(seurat_data, paste(args$output, "_data.h5ad", sep=""))
}

if(args$scope){
    print("Exporting results to SCope compatible loom file")
    io$export_scope_loom(                                                          # we save only counts slot from the RNA assay 
        seurat_data,
        paste(args$output, "_data.loom", sep="")
    )
}