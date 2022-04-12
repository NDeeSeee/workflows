#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(Seurat))
suppressMessages(library(modules))
suppressMessages(library(argparse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(debug <- modules::use(file.path(HERE, "modules/debug.R")))
suppressMessages(filter <- modules::use(file.path(HERE, "modules/filter.R")))
suppressMessages(graphics <- modules::use(file.path(HERE, "modules/graphics.R")))
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(qc <- modules::use(file.path(HERE, "modules/qc.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))


export_all_qc_plots <- function(seurat_data, suffix, args){
    Idents(seurat_data) <- "new.ident"                                                                # safety measure
    selected_features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi")
    selected_labels=c("RNA UMIs", "Genes", "Mitochondrial %", "Novelty score")

    qc_metrics_pca <- qc$qc_metrics_pca(
        seurat_data=seurat_data,
        qc_columns=selected_features,
        qc_labels=selected_labels,
        orq_transform=TRUE
    )

    graphics$pca_plot(
        pca_data=qc_metrics_pca,
        pcs=c(1, 2),
        plot_title=paste(
            paste(
                paste0("PC", c(1, 2)),
                collapse=" and "
            ),
            " from PCA of RNA datasets' QC metrics (", suffix, ")", sep=""
        ),
        legend_title="QC metrics",
        color_by="labels",
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "pca", paste(c(1, 2) ,collapse="_"), "qc_mtrcs", sep="_"),
        pdf=args$pdf
    )
    graphics$pca_plot(
        pca_data=qc_metrics_pca,
        pcs=c(2, 3),
        plot_title=paste(
            paste(
                paste0("PC", c(2, 3)),
                collapse=" and "
            ),
            " from PCA of RNA datasets' QC metrics (", suffix, ")", sep=""
        ),
        legend_title="QC metrics",
        color_by="labels",
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "pca", paste(c(2, 3) ,collapse="_"), "qc_mtrcs", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_bar_plot(
        data=seurat_data@meta.data,
        x_axis="new.ident",
        color_by="new.ident",
        x_label="Identity",
        y_label="Cells",
        legend_title="Identity",
        plot_title=paste("Number of cells per RNA dataset (", suffix, ")", sep=""),
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "cell_count", sep="_"),
        pdf=args$pdf
    )

    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="nCount_RNA",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$rnaminumi,
        x_label="RNA UMIs per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("RNA UMIs per cell density (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "rna_umi_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="nFeature_RNA",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$mingenes,
        x_right_intercept=args$maxgenes,
        x_label="Genes per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Genes per cell density (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "gene_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_point_plot(
        data=seurat_data@meta.data,
        x_axis="nCount_RNA",
        y_axis="nFeature_RNA",
        facet_by="new.ident",
        x_left_intercept=args$rnaminumi,
        y_low_intercept=args$mingenes,
        y_high_intercept=args$maxgenes,
        color_by="mito_percentage",
        gradient_colors=c("lightslateblue", "red", "green"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        x_label="RNA UMIs per cell",
        y_label="Genes per cell",
        legend_title="Mitochondrial %",
        plot_title=paste("Genes vs RNA UMIs per cell correlation (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "gene_rna_umi_corr", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="mito_percentage",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$maxmt,
        x_label="Percentage of transcripts mapped to mitochondrial genes per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("The percentage of transcripts mapped to mitochondrial genes per cell density (", suffix, ")", sep=""),
        zoom_on_intercept=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "mito_perc_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$geom_density_plot(
        data=seurat_data@meta.data,
        x_axis="log10_gene_per_log10_umi",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$minnovelty,
        x_label="log10 Genes / log10 RNA UMIs per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Novelty score per cell density (", suffix, ")", sep=""),
        zoom_on_intercept=TRUE,
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "nvlt_score_dnst", sep="_"),
        pdf=args$pdf
    )
    graphics$vln_plot(
        data=seurat_data,
        features=selected_features,
        labels=selected_labels,
        from_meta=TRUE,
        plot_title=paste("QC metrics per cell densities (", suffix, ")", sep=""),
        legend_title="Identity",
        hide_x_text=TRUE,
        pt_size=0,
        combine_guides="collect",
        palette_colors=graphics$D40_COLORS,
        rootname=paste(args$output, suffix, "qc_mtrcs", sep="_"),
        pdf=args$pdf
    )

    if (seurat_data@meta.data$new.ident != seurat_data@meta.data$condition){
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="nCount_RNA",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=args$rnaminumi,
            x_label="RNA UMIs per cell",
            y_label="Density",
            legend_title="Identity",
            plot_title=paste("Split by grouping condition RNA UMIs per cell density (", suffix, ")", sep=""),
            scale_x_log10=TRUE,
            zoom_on_intercept=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "rna_umi_dnst_spl_by_cond", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="nFeature_RNA",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=args$mingenes,
            x_right_intercept=args$maxgenes,
            x_label="Genes per cell",
            y_label="Density",
            legend_title="Identity",
            plot_title=paste("Split by grouping condition genes per cell density (", suffix, ")", sep=""),
            scale_x_log10=TRUE,
            zoom_on_intercept=TRUE,
            show_ranked=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "gene_dnst_spl_by_cond", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="mito_percentage",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=args$maxmt,
            x_label="Percentage of transcripts mapped to mitochondrial genes per cell",
            y_label="Density",
            legend_title="Identity",
            plot_title=paste("Split by grouping condition the percentage of transcripts mapped to mitochondrial genes per cell density (", suffix, ")", sep=""),
            zoom_on_intercept=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "mito_perc_dnst_spl_by_cond", sep="_"),
            pdf=args$pdf
        )
        graphics$geom_density_plot(
            data=seurat_data@meta.data,
            x_axis="log10_gene_per_log10_umi",
            color_by="new.ident",
            facet_by="condition",
            x_left_intercept=args$minnovelty,
            x_label="log10 Genes / log10 RNA UMIs per cell",
            y_label="Density",
            legend_title="Identity",
            plot_title=paste("Split by grouping condition the novelty score per cell density (", suffix, ")", sep=""),
            zoom_on_intercept=TRUE,
            palette_colors=graphics$D40_COLORS,
            rootname=paste(args$output, suffix, "nvlt_score_dnst_spl_by_cond", sep="_"),
            pdf=args$pdf
        )
    }
}

get_args <- function(){
    parser <- ArgumentParser(description="Seurat RNA Filtering Analysis")
    parser$add_argument(
        "--mex",
        help=paste(
            "Path to the folder with feature-barcode matrix from Cell Ranger Count/Aggregate",
            "experiment in MEX format. If multiple locations provided data is assumed to be not",
            "aggregated (outputs from multiple Cell Ranger Count experiments) and will be merged",
            "before the analysis."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--identity",
        help=paste(
            "Path to the metadata TSV/CSV file to set the datasets identities. If --mex points to",
            "the Cell Ranger Aggregate outputs, the aggregation.csv file can be used. In case of",
            "using feature-barcode matrices from a single or multiple Cell Ranger Count experiments",
            "the file with identities should include at least one column - 'library_id', and a row",
            "with aliases per each experiment from the --mex input. The order of rows should correspond",
            "to the order of feature-barcode matrices provided in the --mex parameter."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--grouping",
        help=paste(
            "Path to the TSV/CSV file to define datasets grouping. First column -",
            "'library_id' with the values provided in the same order as in the",
            "correspondent column of the --identity file, second column 'condition'.",
            "Default: each dataset is assigned to a separate group."
        ),
        type="character"
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
        "--rnamincells",
        help=paste(
            "Include only RNA features detected in at least this many cells. Ignored when",
            "--mex points to the feature-barcode matrices from the multiple Cell Ranger",
            "Count experiments.",
            "Default: 5 (applied to all datasets)"
        ),
        type="integer", default=5
    )
    parser$add_argument(
        "--mingenes",
        help=paste(
            "Include cells where at least this many RNA features are detected.",
            "If multiple values provided, each of them will be applied to the",
            "correspondent dataset from the --mex input based on the --identity",
            "file.",
            "Default: 250 (applied to all datasets)"
        ),
        type="integer", default=250, nargs="*"
    )
    parser$add_argument(
        "--maxgenes",
        help=paste(
            "Include cells with the number of RNA features not bigger than this value.",
            "If multiple values provided, each of them will be applied to the correspondent",
            "dataset from the --mex input based on the --identity file.",
            "Default: 5000 (applied to all datasets)"
        ),
        type="integer", default=5000, nargs="*"
    )
    parser$add_argument(
        "--rnaminumi",
        help=paste(
            "Include cells where at least this many RNA UMIs (transcripts) are detected.",
            "If multiple values provided, each of them will be applied to the correspondent",
            "dataset from the --mex input based on the --identity file.",
            "Default: 500 (applied to all datasets)"
        ),
        type="integer", default=500, nargs="*"
    )
    parser$add_argument(
        "--minnovelty",
        help=paste(
            "Include cells with the novelty score not lower than this value, calculated for",
            "RNA as log10(genes)/log10(UMIs). If multiple values provided, each of them will",
            "be applied to the correspondent dataset from the --mex input based on the",
            "--identity file.",
            "Default: 0.8 (applied to all datasets)"
        ),
        type="double", default=0.8, nargs="*"
    )
    parser$add_argument(
        "--mitopattern",
        help=paste(
            "Regex pattern to identify mitochondrial RNA features.",
            "Default: '^Mt-'"
        ),
        type="character", default="^Mt-"
    )
    parser$add_argument(
        "--maxmt",
        help=paste(
            "Include cells with the percentage of RNA transcripts mapped to mitochondrial",
            "genes not bigger than this value.",
            "Default: 5 (applied to all datasets)"
        ),
        type="double", default=5
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

print(
    paste(
        "Setting parallelization to", args$cpus, "cores, and", args$memory,
        "GB of memory allowed to be shared between the processes"
    )
)
prod$parallel(args)

print(paste("Loading datasets identities from", args$identity))
cell_identity_data <- io$load_cell_identity_data(args$identity)

print(paste("Loading datasets grouping from", args$grouping))
grouping_data <- io$load_grouping_data(args$grouping, cell_identity_data)

print("Loading feature-barcode matrices from:")
for (location in args$mex){print(location)}

seurat_data <- io$load_10x_rna_data(                                                    # identities are set to the "new.ident" column
    args=args,
    cell_identity_data=cell_identity_data,
    grouping_data=grouping_data
)
debug$print_info(seurat_data, args)

print("Adjusting input parameters")
idents_count <- length(unique(as.vector(as.character(Idents(seurat_data)))))
for (key in names(args)){
    if (key %in% c("mingenes", "maxgenes", "rnaminumi", "minnovelty")){
        if (length(args[[key]]) != 1 && length(args[[key]]) != idents_count){
            print(paste("Filtering parameter", key, "has an ambiguous size. Exiting"))
            quit(save="no", status=1, runLast=FALSE)
        }
        if (length(args[[key]]) == 1){
            print(paste("Extending filtering parameter", key, "to have a proper size"))
            args[[key]] <- rep(args[[key]][1], idents_count)
        }
    }
}
print("Adjusted parameters")
print(args)

print(paste("Loading barcodes of interest from", args$barcodes))
barcodes_data <- io$load_barcodes_data(args$barcodes, seurat_data)
print("Applying cell filters based on the loaded barcodes of interest")
seurat_data <- filter$apply_cell_filters(seurat_data, barcodes_data)
debug$print_info(seurat_data, args)

print("Adding RNA QC metrics for not filtered datasets")
seurat_data <- qc$add_rna_qc_metrics(seurat_data, args)
debug$print_info(seurat_data, args)

export_all_qc_plots(
    seurat_data=seurat_data,
    suffix="raw",
    args=args
)

print("Applying filters based on RNA QC metrics")
seurat_data <- filter$apply_rna_qc_filters(seurat_data, cell_identity_data, args)          # cleans up all reductions
debug$print_info(seurat_data, args)

export_all_qc_plots(                                                                       # after all filters have been applied
    seurat_data=seurat_data,
    suffix="fltr",
    args=args
)

DefaultAssay(seurat_data) <- "RNA"                                                         # better to stick to RNA assay by default https://www.biostars.org/p/395951/#395954 
print("Exporting results to RDS file")
io$export_rds(seurat_data, paste(args$output, "_fltr_data.rds", sep=""))
if(args$h5seurat){
    print("Exporting results to h5seurat file")
    io$export_h5seurat(seurat_data, paste(args$output, "_fltr_data.h5seurat", sep=""))
}
