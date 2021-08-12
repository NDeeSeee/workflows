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
suppressMessages(library(ggrepel))
suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(EnhancedVolcano))


set_threads <- function (threads) {
    invisible(capture.output(plan("multiprocess", workers=threads)))
    invisible(capture.output(plan()))
    invisible(capture.output(setDTthreads(threads)))
}


get_filtered_data <- function(seurat_data, args){
    if (!is.null(args$groupby) && !is.null(args$select)){
        print(
            paste(
                "Filtering Seurat data to include only",
                paste(args$select, collapse=" "), "values from",
                args$groupby, "metadata column"
            )
        )
        backup_idents <- Idents(seurat_data)
        Idents(seurat_data) <- args$groupby
        seurat_data <- subset(seurat_data, idents=args$select)
        Idents(seurat_data) <- backup_idents
    }
    return (seurat_data)
}


get_diff_expr_genes <- function(seurat_data, args, assay="RNA", slot="data", min_diff_pct=-Inf){
    backup_assay <- DefaultAssay(seurat_data)
    backup_idents <- Idents(seurat_data)
    DefaultAssay(seurat_data) <- assay
    Idents(seurat_data) <- args$splitby
    diff_expr_markers <- FindMarkers(
        seurat_data,
        slot=slot,
        ident.1=args$first,
        ident.2=args$second,
        logfc.threshold=args$minlogfc,
        min.pct=args$minpct,
        test.use=args$testuse,
        min.diff.pct=min_diff_pct,
        base=2,                      # to make sure we use log2 scale
        verbose=FALSE
    ) %>% rownames_to_column(var="gene")
    Idents(seurat_data) <- backup_idents
    DefaultAssay(seurat_data) <- backup_assay
    return (diff_expr_markers)
}


export_volcano_plot <- function(data, rootname, x_axis, y_axis, x_cutoff, y_cutoff, x_label, y_label, plot_title, plot_subtitle, caption, genes=NULL, label_column="gene", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- EnhancedVolcano(
                data,
                lab=data[,label_column],
                x=x_axis,
                y=y_axis,
                FCcutoff=x_cutoff,
                pCutoff=y_cutoff,
                xlab=x_label,
                ylab=y_label,
                selectLab=genes,
                title=plot_title,
                subtitle=plot_subtitle,
                caption=caption,
                labFace="bold",
                labCol="black",
                # boxedLabels=TRUE,
                # drawConnectors=TRUE,
                # widthConnectors=0.2,
                # endsConnectors="last",
                # typeConnectors="open"
            ) +
            scale_y_log10() +
            theme_gray() +
            NoLegend()

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
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export volcano plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


get_file_type <- function (filename) {
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
            print(paste("Export data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data to", location, sep=" "))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description='Runs differential expression analysis for a subset of cells between two selected groups')
    parser$add_argument(
        "--rds",
        help=paste(
            "Path to the RDS file to load Seurat object from. RDS file can be produced by run_seurat.R script."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Field from the Seurat object metadata to split cells into groups for differential expression analysis.",
            "Default: condition"
        ),
        type="character", default="condition"
    )
    parser$add_argument(
        "--first",
        help=paste(
            "Value from the column set with --splitby to define a first group of cells."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--second",
        help=paste(
            "Value from the column set with --splitby to define a second group of cells."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Field from the Seurat object metadata to group cells for optional subsetting (for example, clustering or predicted cell type field)."
        ),
        type="character"
    )
    parser$add_argument(
        "--select",
        help=paste(
            "Value(s) from the column set with --groupby to optionally subset cells before running differential expression analysis.",
            "Default: do not subset, use all cells."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--genes",
        help="Genes of interest to label on the generated plots. Default: None",
        type="character", nargs="*"
    )
    parser$add_argument(
        "--topn",
        help=paste(
            "Show N genes with the highest and N genes with the lowest log2 fold",
            "change expression values. Ignored when --genes is provided. Default: 10"
        ),
        type="integer", default=10
    )
    parser$add_argument(
        "--minlogfc",
        help=paste(
            "Include only those genes that on average have the absolute value of log2",
            "fold change expression difference not lower than this value. Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--minpct",
        help=paste(
            "Include only those genes that are detected in not lower than this fraction of cells",
            "in either of the two tested groups. Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--maxpvadj",
        help=paste(
            "Include only those genes for which adjusted P-val is not bigger that this value.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--testuse",
        help="Statistical test to use for differential gene expression analysis. Default: wilcox",
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"),
        type="character", default="wilcox"
    )
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./seurat",
        type="character", default="./seurat"
    )
    parser$add_argument(
        "--threads",
        help="Threads. Default: 1",
        type="integer", default=1
    )
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}


args <- get_args()

cat("Step 0: Runtime configuration\n")

print("Used parameters")
print(args)
print(paste("Setting parallelizations threads to", args$threads))
set_threads(args$threads)

cat("\n\nStep 1: Loading data\n")

print(paste("Loading Seurat data from", args$rds))
seurat_data <- readRDS(args$rds)

cat("\n\nStep 2: Running differential expression analysis\n")

filtered_seurat_data <- get_filtered_data(seurat_data, args)
diff_expr_genes <- get_diff_expr_genes(filtered_seurat_data, args) %>%
                   filter(.$p_val_adj<=args$maxpvadj) %>%
                   arrange(desc(avg_log2FC))
topn_diff_expr_genes <- diff_expr_genes %>%
                        filter(row_number() > max(row_number()) - all_of(args$topn) | row_number() <= all_of(args$topn))
topn_diff_expr_genes <- unique(as.vector(as.character(topn_diff_expr_genes[, "gene"])))

highlight_genes <- topn_diff_expr_genes
if ( !is.null(args$genes) ){
    highlight_genes <- args$genes 
}

# print("Exporting differentially expressed genes heatmaps")
# conditions <- c(ident_1, ident_2)
# for (i in 1:length(conditions)){
#     Idents(seurat_data) <- "condition"
#     condition <- conditions[i]
#     cells <- WhichCells(seurat_data, idents=condition)
#     Idents(seurat_data) <- "new.ident"
#     backup_assay <- DefaultAssay(seurat_data)
#     DefaultAssay(seurat_data) <- "RNA"
#     for (i in 1:length(args$resolution)) {
#         resolution <- args$resolution[i]
#         Idents(seurat_data) <- paste("integrated_snn_res", resolution, sep=".")
#         export_expr_heatmap(
#             data=seurat_data,
#             features=diff_expr_markers$gene,
#             cells=cells,
#             plot_title="Log normalized gene expression heatmap of clustered filtered integrated datasets",
#             matrix_slot="data",
#             palette=c("black", "orange"),
#             rootname=paste(args$output, condition, "expr_heatmap_res", resolution, sep="_"),
#             pdf=args$pdf
#         )
#     }
#     Idents(seurat_data) <- "new.ident"
#     DefaultAssay(seurat_data) <- backup_assay
# }

export_volcano_plot(
    data=diff_expr_genes,
    rootname=paste(args$output, "_diff_expr_genes", sep=""),
    x_axis="avg_log2FC",
    y_axis="p_val_adj",
    x_cutoff=args$minlogfc,
    y_cutoff=args$maxpvadj,
    x_label=bquote(~Log[2] ~ "FC"),
    y_label=bquote(~-Log[10] ~ italic(Padj)),
    plot_title="Differentially expressed genes",
    plot_subtitle=paste0(
        args$first, " vs ", args$second, " for cells split by ", args$splitby,
        " (log2FC >= ", args$minlogfc, ", Padj <= ", args$maxpvadj, ")"
    ),
    caption=paste(
        nrow(diff_expr_genes), "genes",
        ifelse(
            (!is.null(args$groupby) && !is.null(args$select)),
            paste("from", paste(args$select, collapse=" "), "subset(s) of", args$groupby, "metadata column"),
            ""
        )
    ),
    genes=highlight_genes
)

print("Exporting filtered differentially expressed genes")
export_data(diff_expr_genes, paste(args$output, "_diff_expr_genes.tsv", sep=""))