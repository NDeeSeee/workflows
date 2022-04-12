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
suppressMessages(filter <- modules::use(file.path(HERE, "modules/filter.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(ucsc <- modules::use(file.path(HERE, "modules/ucsc.R")))


export_all_dimensionality_plots <- function(seurat_data, args) {
    Idents(seurat_data) <- "new.ident"                                                                                         # safety measure
    selected_features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "S.Score", "G2M.Score")
    selected_labels=c("RNA UMIs", "Genes", "Mitochondrial %", "Novelty score", "S score", "G to M score")

    graphics$elbow_plot(
        data=seurat_data,
        ndims=50,
        reduction="pca",
        plot_title="Elbow plot built from PCA of RNA datasets",
        rootname=paste(args$output, "rna_elbow", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_heatmap(
        data=seurat_data,
        dims=1:50,
        plot_title="Genes per cell expression heatmap built from PCA of RNA datasets (sorted by PC scores)",
        x_label="Cells",
        y_label="Genes",
        cells=500,
        nfeatures=30,
        ncol=3,
        height=6500,
        combine_guides="collect",
        rootname=paste(args$output, "rna_pca_heatmap", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_loadings_plot(
        data=seurat_data,
        dims=1:50,
        plot_title="Loadings plot built from PCA of RNA datasets (includes the most variable genes)",
        x_label="PC scores",
        y_label="Genes",
        nfeatures=30,
        ncol=3,
        height=6500,
        combine_guides="collect",
        rootname=paste(args$output, "rna_pca_loadings", sep="_"),
        pdf=args$pdf
    )
    graphics$corr_plot(
        data=seurat_data,
        reduction="pca",
        highlight_dims=args$dimensions,
        qc_columns=selected_features,
        qc_labels=selected_labels,
        plot_title="Correlation plots between QC metrics and principal components from PCA of RNA datasets",
        combine_guides="collect",
        rootname=paste(args$output, "rna_qc_dim_corr", sep="_"),
        pdf=args$pdf
    )
    graphics$feature_plot(
        data=seurat_data,
        features=selected_features,
        labels=selected_labels,
        from_meta=TRUE,
        reduction="rnaumap",
        plot_title="QC metrics on UMAP projected PCA of RNA datasets",
        label=FALSE,
        alpha=0.4,
        max_cutoff="q99",                                                                   # to prevent outlier cells to distort coloring
        combine_guides="keep",
        rootname=paste(args$output, "rna_umap_qc_mtrcs", sep="_"),
        pdf=args$pdf
    )

    graphics$dim_plot(
        data=seurat_data,
        reduction="pca",
        plot_title="PCA of RNA datasets",
        legend_title="Identity",
        group_by="new.ident",
        label=FALSE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, "rna_pca", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title="UMAP projected PCA of RNA datasets",
        legend_title="Identity",
        group_by="new.ident",
        label=FALSE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, "rna_umap", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="pca",
        plot_title="Split by cell cycle phase PCA of RNA datasets",
        legend_title="Identity",
        group_by="new.ident",
        split_by="Phase",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        width=1200,
        height=400,
        rootname=paste(args$output, "rna_pca_spl_by_ph", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title="Split by cell cycle phase UMAP projected PCA of RNA datasets",
        legend_title="Identity",
        group_by="new.ident",
        split_by="Phase",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        width=1200,
        height=400,
        rootname=paste(args$output, "rna_umap_spl_by_ph", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="pca",
        plot_title="Split by quartiles the percentage of transcripts mapped to mitochondrial genes on PCA of RNA datasets",
        legend_title="Identity",
        group_by="new.ident",
        split_by="quartile_mito_percentage",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        width=1200,
        height=400,
        rootname=paste(args$output, "rna_pca_spl_by_mito_qrtl", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title="Split by quartiles the percentage of transcripts mapped to mitochondrial genes on UMAP projected PCA of RNA datasets",
        legend_title="Identity",
        group_by="new.ident",
        split_by="quartile_mito_percentage",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        width=1200,
        height=400,
        rootname=paste(args$output, "rna_umap_spl_by_mito_qrtl", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="pca",
        plot_title="Split by quartiles the UMIs per cell counts on PCA of RNA datasets",
        legend_title="Identity",
        group_by="new.ident",
        split_by="quartile_nCount_RNA",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        width=1200,
        height=400,
        rootname=paste(args$output, "rna_pca_spl_by_umi_qrtl", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title="Split by quartiles the UMIs per cell counts on UMAP projected PCA of RNA datasets",
        legend_title="Identity",
        group_by="new.ident",
        split_by="quartile_nCount_RNA",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        width=1200,
        height=400,
        rootname=paste(args$output, "rna_umap_spl_by_umi_qrtl", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="pca",
        plot_title="Split by quartiles the genes per cell counts on PCA of RNA datasets",
        legend_title="Identity",
        group_by="new.ident",
        split_by="quartile_nFeature_RNA",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        width=1200,
        height=400,
        rootname=paste(args$output, "rna_pca_spl_by_gene_qrtl", sep="_"),
        pdf=args$pdf
    )
    graphics$dim_plot(
        data=seurat_data,
        reduction="rnaumap",
        plot_title="Split by quartiles the genes per cell counts on UMAP projected PCA of RNA datasets",
        legend_title="Identity",
        group_by="new.ident",
        split_by="quartile_nFeature_RNA",
        label=FALSE,
        alpha=0.5,
        palette_colors=graphics$D40_COLORS,
        width=1200,
        height=400,
        rootname=paste(args$output, "rna_umap_spl_by_gene_qrtl", sep="_"),
        pdf=args$pdf
    )

    if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
        graphics$dim_plot(
            data=seurat_data,
            reduction="pca",
            plot_title="Split by identity PCA of RNA datasets",
            legend_title="Identity",
            group_by="new.ident",
            split_by="new.ident",
            label=FALSE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, "rna_pca_spl_by_idnt", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="Split by identity UMAP projected PCA of RNA datasets",
            legend_title="Identity",
            group_by="new.ident",
            split_by="new.ident",
            label=FALSE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, "rna_umap_spl_by_idnt", sep="_"),
            pdf=args$pdf
        )
    }

    if (seurat_data@meta.data$new.ident != seurat_data@meta.data$condition){
        graphics$dim_plot(
            data=seurat_data,
            reduction="pca",
            plot_title="Split by grouping condition PCA of RNA datasets",
            legend_title="Identity",
            group_by="new.ident",
            split_by="condition",
            label=FALSE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, "rna_pca_spl_by_cond", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="Split by grouping condition UMAP projected PCA of RNA datasets",
            legend_title="Identity",
            group_by="new.ident",
            split_by="condition",
            label=FALSE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, "rna_umap_spl_by_cond", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="pca",
            plot_title="Grouped by condition split by cell cycle phase PCA of RNA datasets",
            legend_title="Condition",
            group_by="condition",
            split_by="Phase",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            width=1200,
            height=400,
            rootname=paste(args$output, "rna_pca_gr_by_cond_spl_by_ph", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="Grouped by condition split by cell cycle phase UMAP projected PCA of RNA datasets",
            legend_title="Condition",
            group_by="condition",
            split_by="Phase",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            width=1200,
            height=400,
            rootname=paste(args$output, "rna_umap_gr_by_cond_spl_by_ph", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="pca",
            plot_title="Grouped by condition split by quartiles the percentage of transcripts mapped to mitochondrial genes on PCA of RNA datasets",
            legend_title="Condition",
            group_by="condition",
            split_by="quartile_mito_percentage",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            width=1200,
            height=400,
            rootname=paste(args$output, "rna_pca_gr_by_cond_spl_by_mito_qrtl", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="Grouped by condition split by quartiles the percentage of transcripts mapped to mitochondrial genes on UMAP projected PCA of RNA datasets",
            legend_title="Condition",
            group_by="condition",
            split_by="quartile_mito_percentage",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            width=1200,
            height=400,
            rootname=paste(args$output, "rna_umap_gr_by_cond_spl_by_mito_qrtl", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="pca",
            plot_title="Grouped by condition split by quartiles the UMIs per cell counts on PCA of RNA datasets",
            legend_title="Condition",
            group_by="condition",
            split_by="quartile_nCount_RNA",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            width=1200,
            height=400,
            rootname=paste(args$output, "rna_pca_gr_by_cond_spl_by_umi_qrtl", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="Grouped by condition split by quartiles the UMIs per cell counts on UMAP projected PCA of RNA datasets",
            legend_title="Condition",
            group_by="condition",
            split_by="quartile_nCount_RNA",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            width=1200,
            height=400,
            rootname=paste(args$output, "rna_umap_gr_by_cond_spl_by_umi_qrtl", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="pca",
            plot_title="Grouped by condition split by quartiles the genes per cell counts on PCA of RNA datasets",
            legend_title="Condition",
            group_by="condition",
            split_by="quartile_nFeature_RNA",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            width=1200,
            height=400,
            rootname=paste(args$output, "rna_pca_gr_by_cond_spl_by_gene_qrtl", sep="_"),
            pdf=args$pdf
        )
        graphics$dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title="Grouped by condition split by quartiles the genes per cell counts on UMAP projected PCA of RNA datasets",
            legend_title="Condition",
            group_by="condition",
            split_by="quartile_nFeature_RNA",
            label=FALSE,
            alpha=0.5,
            palette_colors=graphics$D40_COLORS,
            width=1200,
            height=400,
            rootname=paste(args$output, "rna_umap_gr_by_cond_spl_by_gene_qrtl", sep="_"),
            pdf=args$pdf
        )
    } 
}


get_args <- function(){
    parser <- ArgumentParser(description="Seurat RNA Dimensionality Reduction Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load filtered Seurat object from. This file",
            "can be produced by sc_multiome_filter.R, sc_rna_filter.R or sc_atac_reduce.R",
            "scripts. It is mandatory to have gene expression information stored in the RNA assay."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--barcodes",
        help=paste(
            "Path to the headerless TSV/CSV file with the list of barcodes to select",
            "cells of interest (one barcode per line). Prefilters input feature-barcode",
            "matrix to include only selected cells.",
            "Default: use all cells."
        ),
        type="character"
    )
    parser$add_argument(
        "--cellcycle",
        help=paste(
            "Path to the TSV/CSV file with cell cycle data. First column - 'phase',",
            "second column 'gene_id'. If Seurat object loaded from --query parameter",
            "already includes cell cycle scores (S.Score G2M.Score metatada columns),",
            "they will be removed.",
            "Default: skip cell cycle score assignment."
        ),
        type="character"
    )
    parser$add_argument(
        "--rnanorm",
        help=paste(
            "Normalization method to be used for RNA datasets. If provided Seurat object",
            "includes multiple datasets, normalization will be run independently for each",
            "of them, unless integration is disabled with --ntgr set to 'none'",
            "Default: sct"
        ),
        type="character",
        default="sct",
        choices=c("sct", "log", "sctglm")
    )
    parser$add_argument(
        "--ntgr",
        help=paste(
            "Integration method for the joint analysis of multiple RNA datasets present",
            "in the provided with --query Seurat object. Automatically set to 'none' if",
            "only one identity is found.",
            "Default: seurat"
        ),
        type="character",
        default="seurat",
        choices=c("seurat", "none")
    )
    parser$add_argument(
        "--highvarrna",
        help=paste(
            "Number of highly variable RNA features to detect. Used for RNA datasets",
            "integration, scaling, and dimensionality reduction.",
            "Default: 3000"
        ),
        type="integer", default=3000
    )
    parser$add_argument(
        "--regressmt",
        help=paste(
            "Regress mitochondrial genes expression as a confounding source of variation.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--regressrnaumi",
        help=paste(
            "Regress RNA UMIs per cell as a confounding source of variation.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--regressgenes",
        help=paste(
            "Regress genes per cell as a confounding source of variation.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--regresscellcycle",
        help=paste(
            "Regress cell cycle as a confounding source of variation.",
            "Ignored if --cellcycle is not provided.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to use in RNA UMAP projection (from 1 to 50). If single",
            "value N is provided, use from 1 to N PCs. If multiple values are",
            "provided, subset to only selected PCs.",
            "Default: from 1 to 10"
        ),
        type="integer", default=10, nargs="*"
    )
    parser$add_argument(
        "--uspread",
        help=paste(
            "The effective scale of embedded points on UMAP. In combination with",
            "--mindist this determines how clustered/clumped the embedded points are.",
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
        "--lowmem",
        help=paste(
            "Attempts to minimize RAM usage when integrating multiple datasets",
            "with SCTransform algorithm (slows down the computation). Ignored if",
            "--ntgr is not set to 'seurat' or if --rnanorm is not set to either",
            "'sct' or 'sctglm'.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./sc",
        type="character", default="./sc"
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
if (length(args$dimensions) == 1) {
    print("Adjusting --dimensions parameter as only a single value was provided")
    args$dimensions <- c(1:args$dimensions[1])
    print(paste("--dimensions was adjusted to", paste(args$dimensions, collapse=", ")))
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

print(paste("Loading cell cycle data from", args$cellcycle))
cell_cycle_data <- io$load_cell_cycle_data(args$cellcycle)
if(!is.null(cell_cycle_data) && any(c("S.Score", "G2M.Score") %in% colnames(seurat_data@meta.data))){
    print(
        paste(
            "Removing S.Score and G2M.Score columns from",
            "the metadata of the loaded Seurat object",
            "as cell cycle scores will be re-assigned."
        )
    )
    seurat_data[["S.Score"]] <- NULL
    seurat_data[["G2M.Score"]] <- NULL
}

print(paste("Loading barcodes of interest from", args$barcodes))
barcodes_data <- io$load_barcodes_data(args$barcodes, seurat_data)
print("Applying cell filters based on the loaded barcodes of interest")
seurat_data <- filter$apply_cell_filters(seurat_data, barcodes_data)
debug$print_info(seurat_data, args)

print("Running RNA analysis")
seurat_data <- analyses$rna_analyze(seurat_data, args, cell_cycle_data)   # adds "pca" and "rnaumap" reductions
seurat_data <- filter$collapse_fragments_list(seurat_data)                # collapse repetitive fragments if ATAC assay was present in the Seurat object and was splitted
debug$print_info(seurat_data, args)

print("Quantifying QC metrics to evaluate unwanted sources of variation")
seurat_data <- qc$quartile_qc_metrics(
    seurat_data=seurat_data,
    features=c("nCount_RNA", "nFeature_RNA", "mito_percentage"),          # use only those QC metrics that we allow to regress out
    prefix="quartile"                                                     # we will use this prefix in ucsc$export_cellbrowser function
)
debug$print_info(seurat_data, args)

export_all_dimensionality_plots(
    seurat_data=seurat_data,
    args=args
)

if(args$cbbuild){
    print("Exporting RNA assay to UCSC Cellbrowser")
    ucsc$export_cellbrowser(
        seurat_data=seurat_data,
        assay="RNA",
        slot="counts",
        short_label="RNA",
        rootname=paste(args$output, "_cellbrowser", sep=""),
    )
}

DefaultAssay(seurat_data) <- "RNA"                                                         # better to stick to RNA assay by default https://www.biostars.org/p/395951/#395954 
print("Exporting results to RDS file")
io$export_rds(seurat_data, paste(args$output, "_rdcd_data.rds", sep=""))
if(args$h5seurat){
    print("Exporting results to h5seurat file")
    io$export_h5seurat(seurat_data, paste(args$output, "_rdcd_data.h5seurat", sep=""))
}