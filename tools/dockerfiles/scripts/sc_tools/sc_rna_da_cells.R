#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(knitr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(stringr))
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
suppressMessages(logger <- modules::use(file.path(HERE, "modules/logger.R")))

## ----
export_plots <- function(seurat_data, da_cells, pred_thresholds, crnt_thresholds, args) {
    Idents(seurat_data) <- "new.ident"                                                            # safety measure

    datasets_count <- length(unique(as.vector(as.character(seurat_data@meta.data$new.ident))))    # will always be at least 2
    min_dataset_size <- analyses$get_min_ident_size(                                              # all datasets always belong to one of the groups we compare
        SplitObject(seurat_data, split.by="new.ident")
    )
    print(paste("Downsampling to", min_dataset_size, "cells per datasets"))
    downsampled_seurat_data <- subset(seurat_data, downsample=min_dataset_size)                   # downsample per "new.ident"

    graphics$composition_box_plot(
        data=downsampled_seurat_data,                                                             # for a box plot we don't need to downsample per tested condition
        plot_title="Composition box plot colored by tested condition",
        plot_subtitle=paste(
            "Split by cluster;",
            "downsampled to", min_dataset_size,
            "cells per dataset."
        ),
        legend_title="Tested\ncondition",
        group_by=args$splitby,
        split_by=args$groupby,
        stats_by="new.ident",
        x_label="Cluster",
        y_label="Cell counts",
        palette_colors=c(graphics$DOWN_COLOR, graphics$UP_COLOR),
        theme=args$theme,
        rootname=paste(args$output, "cmp_bp_gr_tst_spl_clst", sep="_"),
        pdf=args$pdf
    )

    Idents(downsampled_seurat_data) <- args$splitby
    min_split_size <- analyses$get_min_ident_size(                                                # these are always two groups we compare
        SplitObject(downsampled_seurat_data, split.by=args$splitby)
    )
    print(paste("Downsampling to", min_split_size, "cells per tested condition"))
    downsampled_seurat_data <- subset(downsampled_seurat_data, downsample=min_split_size)         # downsample per tested condition
    Idents(downsampled_seurat_data) <- "new.ident"

    graphics$dim_plot(
        data=downsampled_seurat_data,
        reduction=args$reduction,
        plot_title="UMAP colored by tested condition",
        plot_subtitle=paste(
            "First downsampled to", min_dataset_size,
            "cells per dataset,",
            "then downsampled to", min_split_size,
            "cells per tested condition."
        ),
        legend_title="Tested\ncondition",
        group_by=args$splitby,
        show_density=FALSE,
        label=FALSE,
        pt_size=0.5,
        palette_colors=c(graphics$DOWN_COLOR, graphics$UP_COLOR),
        theme=args$theme,
        rootname=paste(args$output, "umap_gr_tst", sep="_"),
        pdf=args$pdf
    )

    graphics$feature_plot(
        data=seurat_data,
        features="da_score",
        labels=NULL,                                                                                # don't want to show the same text twice
        from_meta=TRUE,
        reduction=args$reduction,
        plot_title="UMAP colored by differential abundance score, categorical scale",
        plot_subtitle=paste0(
            "All cells; ",
            "accepted ranges for differential abundance score: ",
            "[-1, ", round(crnt_thresholds[2], 2), "] and [", round(crnt_thresholds[1], 2), ", 1]"
        ),
        legend_title="Differential\nabundance\nscore",
        label=FALSE,
        pt_size=0.5,
        order=FALSE,
        color_limits=c(-1, 1),
        gradient_colors = c(
            graphics$DOWN_COLOR,
            graphics$DOWN_COLOR,
            graphics$NA_COLOR,
            graphics$NA_COLOR,
            graphics$UP_COLOR,
            graphics$UP_COLOR
        ),
        color_scales = c(
            -1,
            crnt_thresholds[2],
            crnt_thresholds[2]+0.001,
            crnt_thresholds[1]-0.001,
            crnt_thresholds[1],
            1
        ),
        color_breaks = round(c(-1, crnt_thresholds[2], 0, crnt_thresholds[1], 1), 2),
        theme=args$theme,
        rootname=paste(args$output, "umap_da_scr_ctg", sep="_"),
        pdf=args$pdf
    )

    graphics$feature_plot(
        data=seurat_data,
        features="da_score",
        labels=NULL,                                                                    # we already have dataset names on the right side of each plot
        from_meta=TRUE,
        reduction=args$reduction,
        plot_title="UMAP colored by differential abundance score, continuous scale",
        plot_subtitle="All cells",
        legend_title="Differential\nabundance\nscore",
        label=FALSE,
        pt_size=0.5,
        order=FALSE,                                                                    # otherwise white will be on top of blue
        color_limits=c(-1, 1),
        gradient_colors=c(graphics$DOWN_COLOR, graphics$NA_COLOR, graphics$UP_COLOR),
        color_scales=c(-1, 1),
        theme=args$theme,
        rootname=paste(args$output, "umap_da_scr_cnt", sep="_"),
        pdf=args$pdf
    )

    graphics$dim_plot(
        data=downsampled_seurat_data,
        reduction=args$reduction,
        plot_title="UMAP colored by cluster",
        plot_subtitle=paste(
            "Split by tested condition;",
            "first downsampled to", min_dataset_size,
            "cells per dataset,",
            "then downsampled to", min_split_size,
            "cells per tested condition."
        ),
        legend_title="Cluster",
        group_by=args$groupby,
        split_by=args$splitby,
        show_density=TRUE,
        label=FALSE,
        label_color="black",
        pt_size=0.5,
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        # ncol=1,
        width=1200,
        height=600,
        rootname=paste(args$output, "umap_gr_clst_spl_tst", sep="_"),
        pdf=args$pdf
    )

    graphics$composition_plot(
        data=downsampled_seurat_data,
        plot_title="Composition plot colored by cluster",
        plot_subtitle=paste0(
            "Split by tested condition;\n",
            "first downsampled to ", min_dataset_size, " cells per dataset,\n",
            "then downsampled to ", min_split_size, " cells per tested condition."
        ),
        legend_title="Cluster",
        group_by=args$groupby,
        split_by=args$splitby,
        x_label="Tested condition",
        y_label="Cell percentage",
        palette_colors=graphics$D40_COLORS,
        theme=args$theme,
        width=600,                                                              # there are always only two categories
        rootname=paste(args$output, "cmp_gr_clst_spl_tst", sep="_"),
        pdf=args$pdf
    )

    graphics$composition_plot(
        data=downsampled_seurat_data,
        plot_title="Composition plot colored by tested condition",
        plot_subtitle=paste(
            "Split by cluster;",
            "first downsampled to", min_dataset_size,
            "cells per dataset,",
            "then downsampled to", min_split_size,
            "cells per tested condition."
        ),
        legend_title="Tested\ncondition",
        group_by=args$splitby,
        split_by=args$groupby,
        bar_position="dodge",
        x_label="Cluster",
        y_label="Cell counts",
        palette_colors=c(graphics$DOWN_COLOR, graphics$UP_COLOR),
        theme=args$theme,
        rootname=paste(args$output, "cmp_gr_tst_spl_clst", sep="_"),
        pdf=args$pdf
    )

    graphics$daseq_permutations(
        data=da_cells,
        plot_title="Estimated thresholds for differential abundance score",
        plot_subtitle=paste(
            "All cells;",
            args$second, "vs", args$first, "comparison."
        ),
        x_label="Ranked cells",
        y_label="Differential abundance score",
        y_intercepts=round(c(pred_thresholds, crnt_thresholds), 2),
        palette_colors=c("grey", "grey", "black", "black"),
        theme=args$theme,
        width=800,
        height=400,
        rootname=paste(args$output, "rank_da_scr", sep="_"),
        pdf=args$pdf
    )

}

## ----
get_args <- function(){
    parser <- ArgumentParser(description="Single-Cell Differential Abundance Analysis")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load Seurat object from. This file should include",
            "genes expression and/or chromatin accessibility information stored in the RNA",
            "and ATAC assays correspondingly. Both dimensionality reductions selected in",
            "the --reduction and --embeddings parameters should be present in the loaded",
            "Seurat object."
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
        "--embeddings",
        help=paste(
            "Dimensionality reduction to extract embeddings",
            "for differential abundance analysis run with DAseq.",
            "Default: automatically selected based on the",
            "--reduction parameter."
        ),
        type="character"
    )
    parser$add_argument(
        "--dimensions",
        help=paste(
            "Dimensionality to be used when running differential",
            "abundance analysis with DAseq (from 1 to 50).",
            "Default: 10"
        ),
        type="integer", default=10
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to optionally extend Seurat object metadata with",
            "categorical values using samples identities. First column - 'library_id'",
            "should correspond to all unique values from the 'new.ident' column of the",
            "loaded Seurat object. If any of the provided in this file columns are already",
            "present in the Seurat object metadata, they will be overwritten.",
            "Default: no extra metadata is added."
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
            "Default: all cells used, no extra metadata is added."
        ),
        type="character"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from the Seurat object metadata to split cells into two groups",
            "to run --second vs --first differential abundance analysis. May include",
            "columns from the extra metadata added with --metadata or --barcodes",
            "parameters."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--first",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby to define",
            "the first group of cells for differential abundance analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--second",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby to define",
            "the second group of cells for differential abundance analysis."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Column from the Seurat object metadata to group cells",
            "by categories, such as clusters, cell types, etc., when",
            "generating UMAP and composition plots."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--ranges",
        help=paste(
            "Minimum and maximum thresholds to filter out cells with",
            "the low (by absolute values) differential abundance scores.",
            "Default: calculated based on the permutation test."
        ),
        type="double", nargs=2
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
        help="Save raw counts from the RNA and/or ATAC assay(s) to h5ad file(s). Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--loupe",
        help=paste(
            "Save raw counts from the RNA assay to Loupe file.",
            "By enabling this feature you accept the End-User",
            "License Agreement available at https://10xgen.com/EULA.",
            "Default: false"
        ),
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
            "saved. If loaded Seurat object doesn't have RNA assay",
            "this parameter will be ignored. Default: false"
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
    parser$add_argument(
        "--seed",
        help="Seed number for random values. Default: 42",
        type="integer", default=42
    )
    args <- parser$parse_args(str_subset(commandArgs(trailingOnly=TRUE), "\\.R$", negate=TRUE))  # to exclude itself when executed from the sc_report_wrapper.R
    logger$setup(
        file.path(dirname(ifelse(args$output == "", "./", args$output)), "error_report.txt"),
        header="Single-Cell Differential Abundance Analysis (sc_rna_da_cells.R)"
    )
    print(args)
    return (args)
}

## ----
args <- get_args()
prod$parallel(args)

## ----
print(paste("Loading Seurat data from", args$query))
seurat_data <- readRDS(args$query)
debug$print_info(seurat_data, args)

## ----
if (!any(c("RNA", "ATAC") %in% names(seurat_data@assays))){
    logger$info(
        paste(
            "Loaded Seurat object includes neither of the",
            "required assays: RNA and/or ATAC.",
            "Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (args$reduction == "refumap"){
    logger$info(
        paste(
            "Reduction produced be the reference mapping",
            "pipeline (refmap) is not supported. Exiting."
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
if (is.null(args$embeddings)){
    print(
        paste(
            "Attempting to automatically define the",
            "--embeddings parameter on the basis of",
            args$reduction
        )
    )
    args$embeddings <- switch(
        args$reduction,
        "rnaumap" = "pca",
        "atacumap"= "atac_lsi",
        "wnnumap" = "spca"
    )
    if (is.null(args$embeddings)){
        logger$info(
            paste(
                "Failed to automatically define the",
                "--embeddings parameter. Exiting."
            )
        )
        quit(save="no", status=1, runLast=FALSE)
    }
    print(
        paste(
            "Selected", args$embeddings, "embeddings."
        )
    )
}

## ----
if (!(args$embeddings %in% names(seurat_data@reductions))){
    logger$info(
        paste0(
            "Loaded Seurat object doesn't include selected ",
            "embeddings ", args$embeddings, ". Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
print("Adjusting --dimensions parameter")
if (args$embeddings == "atac_lsi"){
    if (!is.null(seurat_data@misc$atac_reduce$first_lsi_removed) && seurat_data@misc$atac_reduce$first_lsi_removed){
        args$dimensions <- c(                                                                    # first LSI component has been already removed
            1:min((args$dimensions - 1), ncol(Embeddings(seurat_data, args$embeddings)))         # use min to make sure it's within the max dimensions
        )
    } else {
        args$dimensions <- c(
            2:min(args$dimensions, ncol(Embeddings(seurat_data, args$embeddings)))               # use min to make sure it's within the max dimensions
        )
    }
} else {
    args$dimensions <- c(
        1:min(args$dimensions, ncol(Embeddings(seurat_data, args$embeddings)))                  # use min to make sure it's within the max dimensions
    )
}
print(paste("--dimensions was adjusted to", paste(args$dimensions, collapse=", ")))

## ----
if(!is.null(args$ranges)){
    args$ranges <- sort(args$ranges, decreasing=TRUE)                                           # always (max, min) thresholds
}

## ----
current_assay <- methods::slot(
    seurat_data@reductions[[args$embeddings]],
    name="assay.used"
)
if (!(current_assay %in% names(seurat_data@assays))){
    logger$info(
        paste(
            "Loaded Seurat object doesn't include required",
            current_assay, "assay. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

print(paste("Setting default assay to", current_assay))
DefaultAssay(seurat_data) <- current_assay                 # shouldn't impact anything because we use cell counts, not the data from the assay
debug$print_info(seurat_data, args)

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
        logger$info(paste("Failed to filter Seurat object with error. Exiting.", e))
        quit(save="no", status=1, runLast=FALSE)
    }
)
debug$print_info(seurat_data, args)

## ----
if (!all(table(seurat_data@meta.data[[args$splitby]]) > 0)){                               # check if we accidentally removed cells we want to compare
    logger$info(
        paste(
            "Not enough cells for comparison. Check --splitby,",
            "--first, --second, and --barcodes parameters. Exiting."
        )
    )
    logger$info(table(seurat_data@meta.data[[args$splitby]]))
    quit(save="no", status=1, runLast=FALSE)
}

## ----
if (!all(args$groupby %in% colnames(seurat_data@meta.data))){
    logger$info(
        paste(
            "Loaded Seurat object can't be grouped by",
            args$groupby, "column. Exiting."
        )
    )
    quit(save="no", status=1, runLast=FALSE)
}

## ----
print(
    paste0(
        "Running ", args$second, " vs ", args$first, " differential abundance ",
        "analysis for datasets split by ", args$splitby, " and grouped by ",
        args$groupby, " using ",
        paste(args$dimensions, collapse=","), " dimensions from the ",
        args$embeddings, " dimensionality reduction.",
        ifelse(
            !is.null(args$ranges),
            paste0(" User provided thresholds for DA scores: ", paste(args$ranges, collapse=",")),
            ""
        )
    )
)
seurat_data@meta.data[[args$splitby]] <- base::factor(                       # to have the proper order on the plots
    seurat_data@meta.data[[args$splitby]],
    levels=c(args$first, args$second)
)
da_results <- analyses$da_analyze(seurat_data, args)                         # will add new metadata column with DA predictions
seurat_data <- da_results$seurat_data                                        # for easy access
da_cells <- da_results$da_cells                                              # we need it only for rand.plot object
pred_thresholds <- da_results$pred_thresholds                                # thresholds identified by DASeq by permutation test
crnt_thresholds <- da_results$crnt_thresholds                                # thresholds used for da_passed_qc column
rm(da_results)                                                               # remove unused data
gc(verbose=FALSE)
debug$print_info(seurat_data, args)

## ----
export_plots(seurat_data, da_cells, pred_thresholds, crnt_thresholds, args)

## ----
if(args$cbbuild){
    if (all(c("RNA", "ATAC") %in% names(seurat_data@assays))){
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            label_field=args$groupby,
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/rna", sep=""),
        )
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            label_field=args$groupby,
            is_nested=TRUE,
            rootname=paste(args$output, "_cellbrowser/atac", sep=""),
        )
    } else if ("RNA" %in% names(seurat_data@assays)){
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="RNA",
            slot="counts",
            short_label="RNA",
            label_field=args$groupby,
            rootname=paste(args$output, "_cellbrowser", sep=""),
        )
    } else {
        ucsc$export_cellbrowser(
            seurat_data=seurat_data,
            assay="ATAC",
            slot="counts",
            short_label="ATAC",
            label_field=args$groupby,
            rootname=paste(args$output, "_cellbrowser", sep="")
        )
    }
}

## ----
io$export_rds(seurat_data, paste(args$output, "_data.rds", sep=""))               # includes only args$first and args$second values from args$splitby column

## ----
if(args$h5seurat){
    io$export_h5seurat(seurat_data, paste(args$output, "_data.h5seurat", sep=""))
}

## ----
if(args$h5ad && ("RNA" %in% names(seurat_data@assays))){
    io$export_h5ad(
        data=seurat_data,
        location=paste(args$output, "_rna_counts.h5ad", sep=""),
        assay="RNA",
        slot="counts"
    )
}

## ----
if(args$h5ad && ("ATAC" %in% names(seurat_data@assays))){
    io$export_h5ad(
        data=seurat_data,
        location=paste(args$output, "_atac_counts.h5ad", sep=""),
        assay="ATAC",
        slot="counts"
    )
}

## ----
if(args$loupe && ("RNA" %in% names(seurat_data@assays))){
    ucsc$export_loupe(
        seurat_data=seurat_data,
        assay="RNA",
        active_cluster=args$groupby,
        rootname=paste0(args$output, "_rna_counts")
    )
}

## ----
if(args$scope && ("RNA" %in% names(seurat_data@assays))){
    io$export_scope_loom(                                                                  # we save only counts slot from the RNA assay 
        seurat_data,
        paste(args$output, "_data.loom", sep="")
    )
}