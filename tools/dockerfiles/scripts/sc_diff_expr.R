#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB should be good by default
options(ggrepel.max.overlaps=Inf)

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
            paste0(
                "Filtering Seurat data to include only '",
                paste(args$select, collapse=" "), "' values from ",
                args$groupby, " metadata column"
            )
        )
        backup_idents <- Idents(seurat_data)
        Idents(seurat_data) <- args$groupby
        seurat_data <- subset(seurat_data, idents=args$select)
        Idents(seurat_data) <- backup_idents
    }
    return (seurat_data)
}


get_avg_expr_data <- function(seurat_data, args, assay="RNA", slot="data"){
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay         # need to explicitely set RNA as scaling is run only on active assay
    avg_expr_data <- AverageExpression(
        seurat_data,
        assays=assay,                          # need only RNA assay
        slot=slot,                             # for slot="data" averaging is done in non-log space
        group.by=args$splitby,
        return.seurat=TRUE,                    # for slot="data" averaged values are saved in "counts", log-normalized - in "data", and scaled - in "scale.data"
        verbose=FALSE
    )
    DefaultAssay(seurat_data) <- backup_assay
    return (avg_expr_data)
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


export_cell_scatter_plot <- function(data, rootname, x_axis, y_axis, plot_title, highlight=NULL, alpha=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            data <- RenameCells(data, new.names=gsub("\\s|\\t|-|\\.|,|\\*|/", "_", Cells(data)))  # otherwise ggplot will fails
            x_axis <- gsub("\\s|\\t|-|\\.|,|\\*|/", "_", x_axis)
            y_axis <- gsub("\\s|\\t|-|\\.|,|\\*|/", "_", y_axis)
            plot <- CellScatter(
                        data,
                        cell1=x_axis,
                        cell2=y_axis,
                        pt.size=2,
                        highlight=highlight
                    ) +
                    ggtitle(plot_title) +
                    theme_gray()

            if (!is.null(highlight)){ plot <- LabelPoints(plot=plot, points=highlight, color="red4", fontface="bold", repel=TRUE) }
            if (!is.null(alpha)) { plot$layers[[1]]$aes_params$alpha <- alpha }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export cell scatter plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export cell scatter plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
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
                        labSize=4,
                        labFace="bold",
                        labCol="red4",
                        colAlpha=0.6,
                        col=c("grey30", "forestgreen", "royalblue", "red"),
                        drawConnectors=TRUE,
                        widthConnectors=0.2
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
        help=paste(
            "Genes of interest to label on the generated plots.",
            "Default: --topn N genes with the highest and the",
            "lowest log2 fold change expression values."
        ),
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

print("Used parameters")
print(args)
print(paste("Setting parallelizations threads to", args$threads))
set_threads(args$threads)

print(paste("Loading Seurat data from", args$rds))
seurat_data <- readRDS(args$rds)

print(
    paste0(
        "Identifying differentially expressed genes between ", args$first, " and ", args$second,
        " for cells split by ", args$splitby, " (log2FC >= ", args$minlogfc, ", Padj <= ", args$maxpvadj, ")"
    )
)
seurat_data <- get_filtered_data(seurat_data, args)
diff_expr_genes <- get_diff_expr_genes(seurat_data, args) %>%
                   filter(.$p_val_adj<=args$maxpvadj) %>%
                   arrange(desc(avg_log2FC))

topn_diff_expr_genes <- diff_expr_genes %>% filter(row_number() > max(row_number()) - all_of(args$topn) | row_number() <= all_of(args$topn))
highlight_genes <- as.vector(as.character(topn_diff_expr_genes[, "gene"]))  # default genes to highlight
if (!is.null(args$genes)){
    print("Check genes of interest to include only those that are differentially expressed")
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(diff_expr_genes[, "gene"]))]
    highlight_genes <- args$genes
}
print(paste("Genes to highlight", paste(highlight_genes, collapse=", ")))

paste("Calculating average gene expression for cells split by", args$splitby)
avg_expr_data <- get_avg_expr_data(seurat_data, args)  # includes only RNA assay. Average normalized expression is saved in "data" slot

export_cell_scatter_plot(                              # uses "data" slot from RNA assay
    data=avg_expr_data,
    rootname=paste(args$output, "avg_gene_expr", sep="_"),
    x_axis=args$first,
    y_axis=args$second,
    highlight=highlight_genes,
    alpha=0.6,
    plot_title=paste("Split by", args$splitby, "log normalized average gene expression"),
    pdf=args$pdf
)

export_volcano_plot(
    data=diff_expr_genes,
    rootname=paste(args$output, "diff_expr_genes", sep="_"),
    x_axis="avg_log2FC",
    y_axis="p_val_adj",
    x_cutoff=args$minlogfc,
    y_cutoff=args$maxpvadj,
    x_label=bquote(~Log[2] ~ "FC"),
    y_label=bquote(~-Log[10] ~ italic(Padj)),
    plot_title="Differentially expressed genes",
    plot_subtitle=paste0(
        args$first, " vs ", args$second, " for cells split by ", args$splitby,
        " (log2FC >= ", args$minlogfc, ", Padj <= ", args$maxpvadj, ")",
        ifelse(
            is.null(args$genes),
            paste0(". Highlight ", 2*args$topn, " genes with the highest abs(log2FC)"),
            ""
        )
    ),
    caption=paste(
        nrow(diff_expr_genes), "genes",
        ifelse(
            (!is.null(args$groupby) && !is.null(args$select)),
            paste0("from '", paste(args$select, collapse=" "), "' subset(s) of ", args$groupby, " metadata column"),
            ""
        )
    ),
    genes=highlight_genes,
    pdf=args$pdf
)

print("Exporting differentially expressed genes")
export_data(diff_expr_genes, paste(args$output, "diff_expr_genes.tsv", sep="_"))