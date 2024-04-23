#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
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


export_all_clustering_plots <- function(seurat_data, args){
    Idents(seurat_data) <- "new.ident"                                                               # safety measure
    datasets_count <- length(unique(as.vector(as.character(seurat_data@meta.data$new.ident))))
    conditions_count <- length(unique(as.vector(as.character(seurat_data@meta.data$condition))))
    not_default_conditions <- all(
        as.vector(as.character(seurat_data@meta.data$new.ident)) != as.vector(as.character(seurat_data@meta.data$condition))
    )

    if (datasets_count > 1){
        Idents(seurat_data) <- "new.ident"
        min_dataset_size <- analyses$get_min_ident_size(SplitObject(seurat_data, split.by="new.ident"))
        print(paste("Downsampling to", min_dataset_size, "cells per datasets"))
        downsampled_per_dataset <- subset(seurat_data, downsample=min_dataset_size)    # downsample per "new.ident"
        if (conditions_count > 1 && not_default_conditions){
            Idents(downsampled_per_dataset) <- "condition"
            min_condition_size <- analyses$get_min_ident_size(SplitObject(downsampled_per_dataset, split.by="condition"))
            print(paste("Additionally downsampling to", min_condition_size, "cells per grouping condition"))
            downsampled_per_condition <- subset(downsampled_per_dataset, downsample=min_condition_size)
            Idents(downsampled_per_dataset) <- "new.ident"
        }
    }

    for (i in 1:length(args$resolution)){
        current_resolution <- args$resolution[i]
        graphics$dim_plot(
            data=seurat_data,
            reduction="atacumap",
            plot_title="UMAP colored by cluster",
            plot_subtitle=paste(
                "All cells;",
                "resolution", current_resolution
            ),
            legend_title="Cluster",
            group_by=paste("atac_res", current_resolution, sep="."),
            label=FALSE,
            label_color="black",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "umap_gr_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        graphics$silhouette_plot(
            data=seurat_data,
            reduction="atac_lsi",
            dims=args$dimensions,
            downsample=500,
            plot_title="Silhouette scores",
            plot_subtitle=paste(
                "All cells;",
                "resolution", current_resolution
            ),
            legend_title="Cluster",
            group_by=paste("atac_res", current_resolution, sep="."),
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "slh_gr_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        if (datasets_count > 1){
            graphics$dim_plot(
                data=downsampled_per_dataset,
                reduction="atacumap",
                plot_title="UMAP colored by cluster",
                plot_subtitle=paste(
                    "Split by dataset;",
                    "downsampled to", min_dataset_size,
                    "cells per dataset;",
                    "resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("atac_res", current_resolution, sep="."),
                split_by="new.ident",
                show_density=TRUE,
                label=FALSE,
                label_color="black",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "umap_gr_clst_spl_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_per_dataset,
                plot_title="Composition plot colored by cluster",
                plot_subtitle=paste(
                    "Split by dataset;",
                    "downsampled to", min_dataset_size,
                    "cells per dataset;",
                    "resolution", current_resolution
                ),
                legend_title="Cluster",
                group_by=paste("atac_res", current_resolution, sep="."),
                split_by="new.ident",
                x_label="Dataset",
                y_label="Cell percentage",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_clst_spl_idnt_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            graphics$composition_plot(
                data=downsampled_per_dataset,
                plot_title="Composition plot colored by dataset",
                plot_subtitle=paste(
                    "Split by cluster;",
                    "downsampled to", min_dataset_size,
                    "cells per dataset;",
                    "resolution", current_resolution
                ),
                legend_title="Dataset",
                group_by="new.ident",
                split_by=paste("atac_res", current_resolution, sep="."),
                bar_position="dodge",
                x_label="Cluster",
                y_label="Cell counts",
                palette_colors=graphics$D40_COLORS,
                theme=args$theme,
                rootname=paste(args$output, "cmp_gr_idnt_spl_clst_res", current_resolution, sep="_"),
                pdf=args$pdf
            )
            if (conditions_count > 1 && not_default_conditions){
                graphics$dim_plot(
                    data=downsampled_per_condition,
                    reduction="atacumap",
                    plot_title="UMAP colored by cluster",
                    plot_subtitle=paste(
                        "Split by grouping condition;",
                        "first downsampled to", min_dataset_size,
                        "cells per dataset,",
                        "then downsampled to", min_condition_size,
                        "cells per grouping condition;",
                        "resolution", current_resolution
                    ),
                    legend_title="Cluster",
                    group_by=paste("atac_res", current_resolution, sep="."),
                    split_by="condition",
                    show_density=TRUE,
                    label=FALSE,
                    label_color="black",
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "umap_gr_clst_spl_cnd_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
                graphics$composition_plot(
                    data=downsampled_per_condition,
                    plot_title="Composition plot colored by cluster",
                    plot_subtitle=paste(
                        "Split by grouping condition;",
                        "first downsampled to", min_dataset_size,
                        "cells per dataset,",
                        "then downsampled to", min_condition_size,
                        "cells per grouping condition;",
                        "resolution", current_resolution
                    ),
                    legend_title="Cluster",
                    group_by=paste("atac_res", current_resolution, sep="."),
                    split_by="condition",
                    x_label="Condition",
                    y_label="Cell percentage",
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "cmp_gr_clst_spl_cnd_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
                graphics$composition_plot(
                    data=downsampled_per_condition,
                    plot_title="Composition plot colored by grouping condition",
                    plot_subtitle=paste(
                        "Split by cluster;",
                        "first downsampled to", min_dataset_size,
                        "cells per dataset,",
                        "then downsampled to", min_condition_size,
                        "cells per grouping condition;",
                        "resolution", current_resolution
                    ),
                    legend_title="Condition",
                    group_by="condition",
                    split_by=paste("atac_res", current_resolution, sep="."),
                    bar_position="dodge",
                    x_label="Cluster",
                    y_label="Cell counts",
                    palette_colors=graphics$D40_COLORS,
                    theme=args$theme,
                    rootname=paste(args$output, "cmp_gr_cnd_spl_clst_res", current_resolution, sep="_"),
                    pdf=args$pdf
                )
            }
        }
    }
    rm(downsampled_per_dataset, downsampled_per_condition)
    gc(verbose=FALSE)
}


export_all_coverage_plots <- function(seurat_data, args) {
    DefaultAssay(seurat_data) <- "ATAC"                                          # safety measure
    Idents(seurat_data) <- "new.ident"                                           # safety measure

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
        for (i in 1:length(args$genes)){
            current_gene <- args$genes[i]
            graphics$coverage_plot(
                data=seurat_data,
                assay="ATAC",
                region=current_gene,
                group_by=paste("atac_res", current_resolution, sep="."),
                plot_title="ATAC fragment coverage",
                plot_subtitle=current_gene,
                idents=NULL,                                                               # to include all values from the default "new.ident" column
                cells=colnames(seurat_data),                                               # limit to only those cells that are in out seurat_data
                features=if("RNA" %in% names(seurat_data@assays)) current_gene else NULL,  # will fail if features are provided without "RNA" assay
                expression_assay="RNA",
                expression_slot="data",                                                    # use scaled counts
                extend_upstream=args$upstream,
                extend_downstream=args$downstream,
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


get_args <- function(){
    parser <- ArgumentParser(description="Single-Cell ATAC-Seq Cluster Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include",
            "chromatin accessibility information stored in the ATAC assay, as well as",
            "'atac_lsi' and 'atacumap' dimensionality reductions applied to that assay."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to use when constructing nearest-neighbor graph before clustering",
            "(from 2 to 50). First LSI component is always excluded unless the provided RDS",
            "file consists of multiple datasets integrated with Harmony.",
            "Default: 10"
        ),
        type="integer", default=10
    )
    parser$add_argument(
        "--ametric",
        help=paste(
            "Distance metric used when constructing nearest-neighbor graph before clustering.",
            "Default: euclidean"
        ),
        type="character", default="euclidean",
        choices=c(
            "euclidean", "cosine", "manhattan", "hamming"
        )
    )
    parser$add_argument(
        "--algorithm",
        help=paste(
            "Algorithm for modularity optimization when running clustering.",
            "Default: slm"
        ),
        type="character", default="slm",
        choices=c(
            "louvain", "mult-louvain", "slm", "leiden"
        )
    )
    parser$add_argument(
        "--resolution",
        help=paste(
            "Clustering resolution applied to the constructed nearest-neighbor graph.",
            "Can be set as an array but only the first item from the list will be used",
            "for cluster labels and peak markers in the UCSC Cell Browser when running",
            "with --cbbuild and --diffpeaks parameters.",
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
            "Genes of interest to build Tn5 insertion frequency plots for the nearest peaks.",
            "If loaded Seurat object includes genes expression information in the RNA assay",
            "it will be additionally shown on the right side of the plots.",
            "Ignored if '--fragments' is not provided.",
            "Default: None"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--upstream",
        help=paste(
            "Number of bases to extend the genome coverage region for",
            "a specific gene upstream. Ignored if --genes or --fragments",
            "parameters are not provided. Default: 2500"
        ),
        type="integer", default=2500
    )
    parser$add_argument(
        "--downstream",
        help=paste(
            "Number of bases to extend the genome coverage region for",
            "a specific gene downstream. Ignored if --genes or --fragments",
            "parameters are not provided. Default: 2500"
        ),
        type="integer", default=2500
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
        "--logfc",
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
        "--minpct",
        help=paste(
            "For differentially accessible peaks identification include only those peaks that",
            "are detected in not lower than this fraction of cells in either of the two tested",
            "clusters. Ignored if '--diffpeaks' is not set.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--testuse",
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
        help="Save raw counts from the ATAC assay to h5ad file. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--cbbuild",
        help="Export results to UCSC Cell Browser. Default: false",
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
    parser$add_argument(
        "--seed",
        help="Seed number for random values. Default: 42",
        type="integer", default=42
    )
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    print(args)
    return (args)
}

args <- get_args()
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

print("Adjusting --dimensions parameter")
if (!is.null(seurat_data@misc$atac_reduce$first_lsi_removed) && seurat_data@misc$atac_reduce$first_lsi_removed){
    args$dimensions <- c(                                                                    # first LSI component has been already removed
        1:min((args$dimensions - 1), ncol(Embeddings(seurat_data, "atac_lsi")))              # use min to make sure it's within the max atac_lsi dimensions
    )
} else {
    args$dimensions <- c(
        2:min(args$dimensions, ncol(Embeddings(seurat_data, "atac_lsi")))                    # use min to make sure it's within the max atac_lsi dimensions
    )
}
print(paste("--dimensions was adjusted to", paste(args$dimensions, collapse=", ")))

if (!is.null(args$fragments)){
    print(paste("Loading ATAC fragments data from", args$fragments))
    seurat_data <- io$replace_fragments(args$fragments, seurat_data)
    debug$print_info(seurat_data, args)
}

print(paste("Clustering ATAC data using", paste(args$dimensions, collapse=", "), "dimensions"))
seurat_data <- analyses$add_clusters(
    seurat_data=seurat_data,
    assay="ATAC",
    graph_name="atac",                          # will be used in all the plot generating functions
    reduction="atac_lsi",
    args=args
)
debug$print_info(seurat_data, args)

export_all_clustering_plots(seurat_data=seurat_data, args=args)

if (!is.null(args$genes)){
    print("Adjusting genes of interest to include only those that are present in the loaded Seurat object")
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(Annotation(seurat_data)$gene_name))]
    if (length(args$genes) == 0){
        print("Neither of the provided genes of interest are present in the loaded Seurat object. Setting back to NULL")
        args$genes <- NULL                                                                                   # all genes were filtered out, so return to NULL
    } else {
        print(paste("--genes was adjusted to", paste(args$genes, collapse=", ")))
    }
}

if (!is.null(args$genes) && !is.null(args$fragments)){
    if ("RNA" %in% names(seurat_data@assays)){
        print("Normalizing counts in RNA assay to show average gene expression alongside the coverage plots")
        DefaultAssay(seurat_data) <- "RNA"
        seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
        DefaultAssay(seurat_data) <- "ATAC"                        # safety measure
    }
    export_all_coverage_plots(
        seurat_data=seurat_data,
        args=args
    )
}

all_markers <- NULL
if (args$diffpeaks){
    print("Identifying differentially accessible peaks between each pair of clusters for all resolutions")
    DefaultAssay(seurat_data) <- "ATAC"                              # safety measure
    all_markers <- analyses$get_markers_by_res(
        seurat_data=seurat_data,
        assay="ATAC",
        resolution_prefix="atac_res",
        latent_vars="nCount_ATAC",                                  # to remove the influence of sequencing depth
        args=args
    )
    if (!is.null(all_markers)){
        io$export_data(
            all_markers,
            paste(args$output, "_peak_markers.tsv", sep="")
        )
    }
}

if(args$cbbuild){
    print("Exporting ATAC assay to UCSC Cellbrowser")
    if(!is.null(all_markers)){
        all_markers <- all_markers %>%
                       dplyr::filter(.$resolution==args$resolution[1]) %>%          # won't fail even if resolution is not present
                       dplyr::select(-c("resolution"))
    }
    print("Reordering reductions to have atacumap on the first place")              # will be shown first in UCSC Cellbrowser
    reduc_names <- names(seurat_data@reductions)
    ordered_reduc_names <- c("atacumap", reduc_names[reduc_names!="atacumap"])      # we checked before that atacumap is present
    seurat_data@reductions <- seurat_data@reductions[ordered_reduc_names]
    debug$print_info(seurat_data, args)
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="ATAC",
        slot="counts",
        short_label="ATAC",
        markers=all_markers,                                                        # can be NULL
        label_field=paste0("Clustering (atac ", args$resolution[1], ")"),           # always use only the first resolution
        palette_colors=graphics$D40_COLORS,                                         # to have colors correspond to the plots
        rootname=paste(args$output, "_cellbrowser", sep="")
    )
}

DefaultAssay(seurat_data) <- "ATAC"
print("Exporting results to RDS file")
io$export_rds(seurat_data, paste(args$output, "_data.rds", sep=""))
if(args$h5seurat){
    print("Exporting results to h5seurat file")
    io$export_h5seurat(seurat_data, paste(args$output, "_data.h5seurat", sep=""))
}

if(args$h5ad){
    print("Exporting ATAC counts to h5ad file")
    io$export_h5ad(
        data=seurat_data,
        location=paste(args$output, "_counts.h5ad", sep=""),
        assay="ATAC",
        slot="counts"
    )
}