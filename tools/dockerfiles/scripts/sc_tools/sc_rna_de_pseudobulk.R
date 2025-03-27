#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(knitr))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(forcats))
suppressMessages(library(stringr))
suppressMessages(library(modules))
suppressMessages(library(argparse))
suppressMessages(library(SummarizedExperiment))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(analyses <- modules::use(file.path(HERE, "modules/analyses.R")))
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(filter <- modules::use(file.path(HERE, "modules/filter.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))
suppressMessages(logger <- modules::use(file.path(HERE, "modules/logger.R")))

## ----
export_raw_plots <- function(seurat_data, args){
    DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    Idents(seurat_data) <- "new.ident"                            # safety measure

    graphics$dim_plot(
        data=seurat_data,
        reduction=args$reduction,
        plot_title="UMAP colored by selected for analysis cells",
        plot_subtitle=paste0(
            "Split by ", args$splitby, "; ",
            ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste(
                    "only cells with",
                    paste(args$subset, collapse=", "),
                    "values in the", args$groupby,
                    "metadata column are selected."
                ),
                "all cells selected."
            )
        ),
        legend_title="Cells",
        split_by=args$splitby,
        group_by=ifelse(                                                          # highlight_group won't work without it
            (!is.null(args$groupby) && !is.null(args$subset)),
            args$groupby,
            args$splitby
        ),
        highlight_group=if(!is.null(args$groupby) && !is.null(args$subset))
                            args$subset
                        else
                            c(args$first, args$second),                           # to highlight all cells
        palette_colors=if(!is.null(args$groupby) && !is.null(args$subset))
                           c(graphics$NA_COLOR, graphics$HIGHLIGHT_COLOR)
                       else
                           c(graphics$HIGHLIGHT_COLOR, graphics$NA_COLOR),        # need to reverse colors
        theme=args$theme,
        rootname=paste(args$output, "umap_spl_tst", sep="_"),
        pdf=args$pdf
    )
    gc(verbose=FALSE)
}

## ----
export_processed_plots <- function(seurat_data, de_results, args){
    DefaultAssay(seurat_data) <- "RNA"                                            # safety measure
    Idents(seurat_data) <- "new.ident"                                            # safety measure

    for (i in 1:length(colnames(de_results$cell$column_metadata))) {
        current_column <- colnames(de_results$cell$column_metadata)[i]
        seurat_data@meta.data[[current_column]] <- factor(
            seurat_data@meta.data[[current_column]],
            levels=unique(                                                        # set levels to the order of appearance in the column_metadata
                as.character(
                    de_results$cell$column_metadata[[current_column]]
                )
            )
        )
    }

    graphics$geom_bar_plot(
        data=seurat_data@meta.data,
        x_axis=ifelse(args$test %in% c("deseq", "lrt"), "new.ident", args$splitby),
        color_by=args$splitby,
        x_label=ifelse(args$test %in% c("deseq", "lrt"), "Dataset", args$splitby),
        y_label="Cell counts",
        legend_title=args$splitby,
        plot_title=ifelse(
            args$test %in% c("deseq", "lrt"),
            "Number of cells per dataset",
            paste(
                "Number of cells per tested",
                "condition defined by", args$splitby,
                "metadata column"
            )
        ),
        plot_subtitle=paste0(
            "Colored by ", args$splitby, "; ",
            ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste(
                    "only cells with",
                    paste(args$subset, collapse=", "),
                    "values in the", args$groupby,
                    "metadata column are selected."
                ),
                "all cells selected."
            )
        ),
        palette_colors=c(graphics$TRUE_COLOR, graphics$FALSE_COLOR),
        theme=args$theme,
        rootname=paste(args$output, "cell_cnts", sep="_"),
        pdf=args$pdf
    )

    if (!is.null(de_results$bulk)){
        graphics$mds_html_plot(
            norm_counts_data=de_results$bulk$counts_data,               # not filtered pseudobulk aggregated normalized counts data in a form of SummarizedExperiment
            rootname=paste(args$output, "mds_plot", sep="_")
        )

        pca_data <- qc$counts_pca(                                      # adds group column to identify the datasets
            SummarizedExperiment::assay(de_results$bulk$counts_data)
        )

        graphics$pca_plot(
            pca_data=pca_data,
            pcs=c(1, 2),
            plot_title="Gene expression PCA",
            plot_subtitle="PC1/PC2",
            legend_title="Dataset",
            color_by="group",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "pca_1_2", sep="_"),
            pdf=args$pdf
        )

        graphics$pca_plot(
            pca_data=pca_data,
            pcs=c(2, 3),
            plot_title="Gene expression PCA",
            plot_subtitle="PC2/PC3",
            legend_title="Dataset",
            color_by="group",
            palette_colors=graphics$D40_COLORS,
            theme=args$theme,
            rootname=paste(args$output, "pca_2_3", sep="_"),
            pdf=args$pdf
        )
    }

    genes_to_highlight <- get_genes_to_highlight(de_results, args)  # may return empty vector
    volcano_plot_data <- de_results$de_genes %>%
                         dplyr::mutate(
                             size_by=dplyr::if_else(
                                 .$padj <= args$padj & abs(.$log2FoldChange) >= args$logfc,
                                 dplyr::if_else(
                                     .$log2FoldChange > 0,
                                     .[[base::paste("pct", args$second, sep="_")]],
                                     .[[base::paste("pct", args$first, sep="_")]]
                                 ),
                                 0
                             ),
                             overlay_by=dplyr::if_else(
                                 .$padj <= args$padj & abs(.$log2FoldChange) >= args$logfc,
                                 dplyr::if_else(
                                     .$log2FoldChange > 0,
                                     .[[base::paste("pct", args$first, sep="_")]],
                                     .[[base::paste("pct", args$second, sep="_")]]
                                 ),
                                 0
                             ),
                             alpha_by=dplyr::if_else(
                                 .$padj <= args$padj & abs(.$log2FoldChange) >= args$logfc,
                                 1,
                                 0.5
                             )
                         ) %>%
                         dplyr::mutate(
                             color_by=dplyr::if_else(
                                 .$padj <= args$padj & abs(.$log2FoldChange) >= args$logfc,
                                 dplyr::if_else(
                                    .$overlay_by <= .$size_by,
                                    dplyr::if_else(
                                        .$log2FoldChange > 0,
                                        "firebrick2",
                                        "dodgerblue"
                                    ),
                                    "#865E96"
                                 ),
                                 graphics$NA_COLOR
                             ),
                             overlay_color_by=dplyr::if_else(
                                 .$padj <= args$padj & abs(.$log2FoldChange) >= args$logfc,
                                 dplyr::if_else(
                                    .$overlay_by > .$size_by,
                                    dplyr::if_else(
                                        .$log2FoldChange > 0,
                                        "dodgerblue",
                                        "firebrick2"
                                    ),
                                    "#865E96"
                                 ),
                                 graphics$NA_COLOR
                             )
                         )

    graphics$volcano_plot(
        data=volcano_plot_data,                                   # this is not filtered differentially expressed features
        x_axis="log2FoldChange",
        y_axis="padj",
        x_cutoff=args$logfc,
        y_cutoff=args$padj,
        x_label="log2 FC",
        y_label="-log10 Padj",
        plot_title="Differentially expressed genes",
        plot_subtitle=paste0(
            "Split by ", args$splitby, " (", args$second, " vs ", args$first, "); ",
            ifelse(
                (!is.null(args$groupby) && !is.null(args$subset)),
                paste(
                    "only cells with", paste(args$subset, collapse=", "),
                    "values in", args$groupby, "metadata column are selected; "
                ),
                "all cells selected; "
            ),
            "use ", args$test, " test; ",
            ifelse(
                (!is.null(args$batchby)),
                paste0("model batch effect by ", args$batchby, "; "),
                ""
            ),
            "adjusted P-value threshold ", args$padj, "; ",
            "log2FC threshold ", args$logfc, "; ",
            ifelse(
                (!is.null(args$genes) && length(args$genes) > 0),
                "label user provided genes",
                paste("label top", length(genes_to_highlight), "genes")
            )
        ),
        caption=paste(nrow(volcano_plot_data), "genes"),
        features=genes_to_highlight,
        size_by="size_by",
        size_from=c(0, 1),
        color_by="color_by",
        alpha_by="alpha_by",
        overlay_by="overlay_by",
        overlay_from=c(0, 1),
        overlay_color_by="overlay_color_by",
        theme=args$theme,
        rootname=paste(args$output, "dxpr_vlcn", sep="_"),
        pdf=args$pdf
    )

    if (!is.null(genes_to_highlight) && length(genes_to_highlight) > 0){
        graphics$dot_plot(
            data=seurat_data,
            features=genes_to_highlight,
            plot_title="Average gene expression",
            plot_subtitle=paste(
                "Split by",
                ifelse(
                    args$test %in% c("deseq", "lrt"),
                    "dataset",
                    args$splitby
                )
            ),
            x_label="Genes",
            y_label=ifelse(args$test %in% c("deseq", "lrt"), "Dataset", args$splitby),
            cluster_idents=FALSE,
            group_by=ifelse(args$test %in% c("deseq", "lrt"), "new.ident", args$splitby),
            scale=args$test %in% c("deseq", "lrt"),
            theme=args$theme,
            rootname=paste(args$output, "xpr_avg", sep="_"),
            pdf=args$pdf
        )
        graphics$vln_plot(
            data=seurat_data,
            features=genes_to_highlight,
            labels=genes_to_highlight,
            plot_title="Gene expression density",
            plot_subtitle=paste("Split by", args$splitby),
            legend_title=args$splitby,
            rotate_labels=TRUE,
            pt_size=0,
            hide_x_text=TRUE,
            group_by=args$splitby,
            combine_guides="collect",
            palette_colors=c(graphics$TRUE_COLOR, graphics$FALSE_COLOR),   # we always split by args$splitby so we have only two categories
            theme=args$theme,
            rootname=paste(args$output, "xpr_dnst", sep="_"),
            pdf=args$pdf
        )
        for (current_gene in genes_to_highlight) {
            graphics$feature_plot(
                data=seurat_data,
                features=current_gene,
                labels=current_gene,
                reduction=args$reduction,
                plot_title="UMAP colored by gene expression",
                plot_subtitle=paste0(
                    "Split by ", args$splitby, "; ",
                    ifelse(
                        (!is.null(args$groupby) && !is.null(args$subset)),
                        paste(
                            "only cells with",
                            paste(args$subset, collapse=", "),
                            "values in the", args$groupby,
                            "metadata column are selected."
                        ),
                        "all cells selected."
                    )
                ),
                legend_title="Expression",
                label=FALSE,
                order=TRUE,
                split_by=args$splitby,
                pt_size=1,                     # need to have it fixed, otherwise dots are different on the splitted plot
                combine_guides="collect",
                theme=args$theme,
                width=800,
                height=400,
                rootname=paste(args$output, "xpr_per_cell", current_gene, sep="_"),
                pdf=args$pdf
            )
        }
    }

    if ("HCL" %in% colnames(de_results$cell$row_metadata)){       # safe to take from cells, as if present in bulk it will be still the same
        split_rows <- forcats::fct_inorder(as.character(de_results$cell$row_metadata[["HCL"]]))
    } else if ("HCL.1" %in% colnames(de_results$cell$row_metadata)){
        split_rows <- forcats::fct_inorder(as.character(de_results$cell$row_metadata[["HCL.1"]]))
    } else {
        split_rows <- NULL
    }

    graphics$feature_heatmap(                                      # we build it only for cells level
        data=seurat_data,
        assay="RNA",                                               # for now we will always use RNA even if SCT may be present
        slot="data",
        cells=rownames(de_results$cell$column_metadata),           # to make sure we use the same number and order of cells as in GCT file
        features=rownames(de_results$cell$row_metadata),           # to make sure we use the same number and order of genes as in GCT file
        split_rows=split_rows,
        show_rownames=FALSE,
        scale_to_max=FALSE,
        scale="row",                                               # will calculate z-score
        heatmap_colors=c("darkblue", "black", "yellow"),
        group_by=colnames(de_results$cell$column_metadata),        # doesn't change the order, only used in annotations
        order_by=NULL,                                             # we don't want to reorder cells by columns defined in the group_by
        palette_colors=graphics$D40_COLORS,
        plot_title=paste(
            "Gene expression heatmap filtered",
            "by P adjusted <=", args$padj, "and",
            "|log2FC| >=", args$logfc,
            switch(
                ifelse(is.null(args$cluster), "none", args$cluster),
                "row" = "(clustered by gene)",
                "column" = "(clustered by dataset)",
                "both" = "(clustered both by gene and dataset)",
                "none" = "(not clustered)"
            )
        ),
        rootname=paste(args$output, "xpr_htmp", sep="_"),
        pdf=args$pdf
    )

    io$export_data(
        de_results$cell$row_metadata %>% tibble::rownames_to_column(var="feature"),
        paste(args$output, "xpr_htmp.tsv", sep="_")
    )

}

## ----
get_genes_to_highlight <- function(de_results, args){
    if (!is.null(args$genes) && length(args$genes) > 0){
        print("Using user provided genes to highlight regardless of the significance score")
        genes_to_highlight <- args$genes[args$genes %in% as.vector(as.character(de_results$de_genes[, "gene"]))]
    } else {
        print(
            paste(
                "Identifying up to 20 the most",
                "differentially expressed genes",
                "with P adjusted <=", args$padj,
                "and |log2FC| >=", args$logfc
            )
        )
        top_up_genes <- de_results$de_genes %>%
                        filter(.$padj <= args$padj) %>%
                        filter(.$log2FoldChange >= args$logfc) %>%
                        arrange(desc(log2FoldChange)) %>%
                        slice_head(n=10)
        top_down_genes <- de_results$de_genes %>%
                          filter(.$padj <= args$padj) %>%
                          filter(.$log2FoldChange <= -args$logfc) %>%
                          arrange(log2FoldChange) %>%
                          slice_head(n=10)
        genes_to_highlight <- c(
            as.vector(as.character(top_up_genes$gene)),
            as.vector(as.character(top_down_genes$gene))
        )
    }
    print(paste("Genes to highlight", paste(genes_to_highlight, collapse=", ")))
    return (genes_to_highlight)
}

## ----
get_args <- function(){
    parser <- ArgumentParser(
        description="Single-Cell RNA-Seq Differential Expression Analysis"
    )
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from.",
            "This file should include genes expression information",
            "stored in the RNA assay. The dimensionality reductions",
            "selected in the --reduction parameter should be present",
            "in the loaded Seurat object."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--reduction",
        help=paste(
            "Dimensionality reduction to be used for generating UMAP plots."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to optionally extend Seurat object metadata with",
            "categorical values using samples identities. First column - library_id",
            "should correspond to all unique values from the new.ident column of the",
            "loaded Seurat object. If any of the provided in this file columns are already",
            "present in the Seurat object metadata, they will be overwritten. When combined",
            "with --barcodes parameter, first the metadata will be extended, then barcode",
            "filtering will be applied. Default: no extra metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--barcodes",
        help=paste(
            "Path to the TSV/CSV file to optionally prefilter and extend Seurat object",
            "metadata by selected barcodes. First column should be named as barcode.",
            "If file includes any other columns they will be added to the Seurat object",
            "metadata ovewriting the existing ones if those are present.",
            "Default: all cells used, no extra metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Column from the Seurat object metadata to group cells for optional subsetting",
            "when combined with --subset parameter. May be one of the extra metadata columns",
            "added with --metadata or --barcodes parameters. Ignored if --subset is not set.",
            "Default: do not subset, include all cells into analysis."
        ),
        type="character"
    )
    parser$add_argument(
        "--subset",
        help=paste(
            "Values from the column set with --groupby parameter to subset cells",
            "before running differential expression analysis. Ignored if --groupby",
            "is not provided. Default: do not subset cells, include all of them."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from the Seurat object metadata to split cells into two groups to",
            "run --second vs --first differential expression analysis. If --test",
            "parameter is set to deseq or deseq-lrt, the --splitby shouldn't put cells",
            "from the same dataset into the different comparison groups. May be one of",
            "the extra metadata columns added with --metadata or --barcodes parameters."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--first",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby parameter",
            "to define the first group of cells for differential expression analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--second",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby parameter",
            "to define the second group of cells for differential expression analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--test",
        help=paste(
            "Test type to use in differential expression analysis. If set to deseq or",
            "deseq-lrt, gene expression will be aggregated to the pseudobulk level per",
            "dataset. For deseq, the pair-wise Wald test will be used. For deseq-lrt, the",
            "reduced formula will look like ~1 if --batchby parameter is omitted or will",
            "be set to ~batchby to exclude the criteria if interest (defined by --splitby).",
            "For all other values of the --test parameter the FindMarkers function will",
            "be used (genes will be prefiltered by minimum percentage >= 0.1 and by minimum",
            "log2FoldChange >= 0.25 before running differential expression analysis).",
            "Default: use FindMarkers with Wilcoxon Rank Sum test."
        ),
        type="character", default="wilcoxon",
        choices=c(
            "wilcoxon",                   # (wilcox) Wilcoxon Rank Sum test
            "likelihood-ratio",           # (bimod) Likelihood-ratio test
            "t-test",                     # (t) Student's t-test
            "negative-binomial",          # (negbinom) Negative Binomial Generalized Linear Model (supports --batchby)
            "poisson",                    # (poisson) Poisson Generalized Linear Model (supports --batchby)
            "logistic-regression",        # (LR) Logistic Regression (supports --batchby)
            "mast",                       # (MAST) MAST package (supports --batchby)
            "deseq",                      # DESeq2 Wald test on pseudobulk aggregated gene expression
            "deseq-lrt"                   # DESeq2 LRT test on pseudobulk aggregated gene expression
        )
    )
    parser$add_argument(
        "--batchby",
        help=paste(
            "Column from the Seurat object metadata to group cells into batches.",
            "If --test is set to deseq or deseq-lrt the --batchby parameter will",
            "be used in the design formula in the following way ~splitby+batchby.",
            "Additionally, the --batchby shouldn't put cells from the same dataset",
            "into the different batches. If --test is set to negative-binomial,",
            "poisson, logistic-regression, or mast it will be used as a latent",
            "variable in the FindMarkers function. Not supported for --test values",
            "equal to wilcoxon, likelihood-ratio, or t-test. May be one of the extra",
            "metadata columns added with --metadata or --barcodes parameters.",
            "Default: do not model batch effect."
        ),
        type="character"
    )
    parser$add_argument(
        "--padj",
        help=paste(
            "In the exploratory visualization part of the analysis output only",
            "differentially expressed genes with adjusted P-value not bigger",
            "than this value. Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--minpct",
        help=paste(
            "Include only those genes that are detected in not lower than this",
            "fraction of cells in either of the two tested conditions.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--logfc",
        help=paste(
            "In the exploratory visualization part of the analysis output only",
            "differentially expressed genes with log2 Fold Change not smaller than",
            "this value. Default: 0.585 (1.5 folds)"
        ),
        type="double", default=0.585
    )
    parser$add_argument(
        "--genes",
        help=paste(
            "Genes of interest to label on the generated plots. Default: top 10",
            "genes with the highest and the lowest log2FoldChange values."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--exclude",
        help=paste(
            "Regex pattern to identify and exclude specific genes from the",
            "differential expression analysis (not case-sensitive). If any",
            "of such genes are provided in the --genes parameter, they will",
            "be excluded from there as well. Default: use all genes"
        ),
        type="character"
    )
    parser$add_argument(
        "--cluster",
        help=paste(
            "Hopach clustering method to be run on the normalized read counts for",
            "the exploratory visualization part of the analysis. Clustering by",
            "column is supported only when --test is set to deseq or deseq-lrt.",
            "Default: do not run clustering"
        ),
        type="character",
        choices=c("row", "column", "both")
    )
    parser$add_argument(
        "--rowdist",
        help=paste(
            "Distance metric for HOPACH row clustering. Ignored if --cluster is set",
            "to column or not provided. Default: cosangle"
        ),
        type="character", default="cosangle",
        choices=c("cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor")
    )
    parser$add_argument(
        "--columndist",
        help=paste(
            "Distance metric for HOPACH column clustering. Ignored if --cluster is set",
            "to row or not provided. Default: euclid"
        ),
        type="character", default="euclid",
        choices=c("cosangle", "abscosangle", "euclid", "abseuclid", "cor", "abscor")
    )
    parser$add_argument(
        "--center",
        help=paste(
            "Apply mean centering for gene expression prior to running",
            "clustering by row. Ignored if --cluster is set to column or",
            "not provided. Default: do not centered"
        ),
        action="store_true"
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
    args <- parser$parse_args(str_subset(commandArgs(trailingOnly=TRUE), "\\.R$", negate=TRUE))  # to exclude itself when executed from the sc_report_wrapper.R
    logger$setup(
        file.path(dirname(ifelse(args$output == "", "./", args$output)), "error_report.txt"),
        header="Single-Cell RNA-Seq Differential Expression Analysis (sc_rna_de_pseudobulk.R)"
    )
    print(args)
    return (args)
}

## ----
args <- get_args()
print("Adjusting --test parameter")
args$test <- switch(
    args$test,
    "wilcoxon"            = "wilcox",
    "likelihood-ratio"    = "bimod",
    "t-test"              = "t",
    "negative-binomial"   = "negbinom",
    "poisson"             = "poisson",
    "logistic-regression" = "LR",
    "mast"                = "MAST",
    "deseq"               = "deseq",
    "deseq-lrt"           = "lrt"
)
print(paste("--test was adjusted to", args$test))
prod$parallel(args)

## ----
print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)
debug$print_info(seurat_data, args)

## ----
if (!("RNA" %in% names(seurat_data@assays))){
    logger$info(
        paste(
            "Loaded Seurat object doesn't include",
            "the required RNA assay. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (!(args$reduction %in% names(seurat_data@reductions))){
    logger$info(
        paste0(
            "Loaded Seurat object doesn't include selected ",
            "reduction ", args$reduction, ". Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if ( !is.null(args$cluster) && (args$cluster %in% c("column", "both")) && !(args$test %in% c("deseq", "lrt")) ){
    logger$info(
        paste(
            "Clustering by column is not supported if gene expression data",
            "are not aggregated to pseudobulk form. Check --test parameter.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (!is.null(args$batchby) && (args$test %in% c("wilcox", "bimod", "t"))){
    logger$info(
        paste(
            "Modeling batch effect is not supported if --test parameter",
            "was set to wilcoxon, likelihood-ratio, or t-test. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
print("Setting default assay to RNA")
DefaultAssay(seurat_data) <- "RNA"

## ----
excluded_genes <- c()
if (!is.null(args$exclude)){
    excluded_genes <- grep(
        args$exclude,
        as.vector(as.character(rownames(seurat_data))),    # with RNA assay set as default the rownames should be genes
        value=TRUE,
        ignore.case=TRUE
    )
    print(
        paste(
            "Based on the pattern", args$exclude, "the following genes will",
            "be excluded from the differential expression analysis:",
            paste(excluded_genes, collapse=", ")
        )
    )
}

## ----
if (!is.null(args$genes) && length(args$genes) > 0){
    print(
        paste(
            "Adjusting genes of interest to include only those genes that are",
            "present in the loaded Seurat object and have not been excluded."
        )
    )
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(rownames(seurat_data)))]     # with RNA assay set as default the rownames should be genes
    args$genes <- args$genes[!(args$genes %in% excluded_genes)]                                  # excluded_genes can be empty list
    print(args$genes)
}

## ----
if (!is.null(args$metadata)){
    print("Extending Seurat object with the extra metadata fields")
    seurat_data <- io$extend_metadata(
        seurat_data=seurat_data,
        location=args$metadata,
        seurat_ref_column="new.ident",
        meta_ref_column="library_id"
    )
    debug$print_info(seurat_data, args)
}

## ----
if (!is.null(args$barcodes)){
    print("Applying cell filters based on the barcodes of interest")
    seurat_data <- io$extend_metadata_by_barcode(seurat_data, args$barcodes, TRUE)         # sets identities to new.ident
    debug$print_info(seurat_data, args)
}

## ----
if (args$test %in% c("deseq", "lrt")){
    # need to make sure that neither --splitby nor --batchby
    # don't put cells from the same dataset into the different
    # groups, because when we use deseq or lrt test we aggregate
    # everything to the pseudobulk form per dataset.

    for (criteria in c(args$splitby, args$batchby)){
        if (!is.null(criteria)){                                                           # because --batchby can be NULL
            cells_counts <- table(
                seurat_data@meta.data$new.ident,
                seurat_data@meta.data[[criteria]]
            )
            if (any(rowSums(cells_counts > 0) > 1)){
                logger$info(
                    paste(
                        "Dividing cells by", criteria, "puts cells",
                        "from the same dataset into the different",
                        "comparison groups or batches, which is not",
                        "supported when --test parameter is set to",
                        "deseq or deseq-lrt. Exiting."
                    )
                )
                logger$info(cells_counts)
                quit(save="no", status=1, runLast=FALSE)
            }
        }
    }
}

## ----
tryCatch(
    expr = {
        print("Subsetting Seurat object to include only cells from the tested conditions")
        seurat_data <- io$apply_metadata_filters(
            seurat_data,
            args$splitby,
            c(args$first, args$second)
        )
    },
    error = function(e){
        logger$info(
            paste(
                "Failed to filter Seurat object by",
                paste(c(args$first, args$second), collapse=", "),
                "values from the", args$splitby, "column.",
                "Exiting.", e
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
)
debug$print_info(seurat_data, args)

## ----
export_raw_plots(seurat_data, args)                                                        # will include only cells from --first and --second

## ----
if(!is.null(args$groupby) && !is.null(args$subset)){
    tryCatch(
        expr = {
            print("Subsetting Seurat object to include only selected groups of cells")
            seurat_data <- io$apply_metadata_filters(
                seurat_data,
                args$groupby,
                args$subset
            )
        },
        error = function(e){
            logger$info(
                paste(
                    "Failed to filter Seurat object by",
                    paste(args$subset, collapse=", "),
                    "values from the", args$subset, "column.",
                    "Exiting.", e
                )
            )
            quit(save="no", status=1, runLast=FALSE)
        }
    )
    debug$print_info(seurat_data, args)
}

## ----
if (!all(table(seurat_data@meta.data[[args$splitby]]) > 0)){                               # check if we accidentally removed cells we want to compare
    logger$info(
        paste(
            "Not enough cells for comparison. Check --groupby",
            "and --subset parameters. Exiting."
        )
    )
    print(table(seurat_data@meta.data[[args$splitby]]))
    quit(save="no", status=1, runLast=FALSE)
}

## ----
print("Normalizing RNA counts after all filters applied")
seurat_data <- NormalizeData(seurat_data, verbose=FALSE)                                   # just in case rerun normalize as we use data slot in FindMarkers function

## ----
de_results <- analyses$rna_de_analyze(seurat_data, args, excluded_genes)
print(head(de_results$de_genes))

## ----
export_processed_plots(seurat_data, de_results, args)

## ----
io$export_data(
    de_results$de_genes,                                                                    # not filtered na-removed differentially expressed genes
    paste0(args$output, "_de_genes.tsv")
)

## ----
if (!is.null(de_results$bulk)){                                                             # looks like we run pseudobulk analysis
    io$export_gct(                                                                          # will be used by GSEA
        counts_mat=SummarizedExperiment::assay(de_results$bulk$counts_data),
        row_metadata=NULL,                                                                  # should be NULL, because de_results$row_metadata is filtered by padj
        col_metadata=de_results$bulk$column_metadata,                                       # includes samples as row names
        location=paste0(args$output, "_bulk_counts.gct")
    )
    io$export_cls(
        categories=de_results$bulk$column_metadata[, args$splitby],                         # to have the same samples order as in GCT file
        paste0(args$output, "_bulk_phntps.cls")
    )
}

## ----
print("Exporting differentially expressed genes as morpheus heatmap")
expression_mat <- t(scale(t(de_results$cell$counts_mat)))                                   # will calculate z-score
expression_limits <- stats::quantile(                                                       # to exclude outliers
    abs(expression_mat), 0.99, na.rm=TRUE, names=FALSE
)

## ----
io$export_gct(                                                                              # will be used by Morpheus
    counts_mat=expression_mat,
    row_metadata=de_results$cell$row_metadata,                                              # includes genes as row names
    col_metadata=de_results$cell$column_metadata,                                           # includes cells as row names
    location=paste0(args$output, "_xpr_htmp.gct")
)

## ----
graphics$morpheus_html_heatmap(
    gct_location=paste0(args$output, "_xpr_htmp.gct"),
    rootname=paste0(args$output, "_xpr_htmp"),
    color_scheme=list(
        scalingMode="fixed",
        stepped=FALSE,
        values=as.list(c(-expression_limits, 0, expression_limits)),
        colors=c("darkblue", "black", "yellow")
    )
)
