#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB should be good by default
options(ggrepel.max.overlaps=Inf)

suppressMessages(library(scran))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(DESeq2))
suppressMessages(library(tibble))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(argparse))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(BiocParallel))
suppressMessages(library(EnhancedVolcano))


set_threads <- function (threads) {
    invisible(capture.output(plan("multiprocess", workers=threads)))
    invisible(capture.output(plan()))
    invisible(capture.output(setDTthreads(threads)))
    register(MulticoreParam(threads))                  # for DESeq2
}


extend_metadata <- function(seurat_data, args) {
    if (!is.null(location)){
        extra_metadata <- read.table(
            args$condition,
            sep=get_file_type(args$condition),
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )
        print(extra_metadata)
        Idents(seurat_data) <- "new.ident"
        identity_data <- unique(as.vector(as.character(Idents(seurat_data))))
        if ( (nrow(extra_metadata) == length(identity_data)) && all(is.element(identity_data, extra_metadata$library_id)) ){
            print(paste("Extra metadata is successfully loaded from ", args$condition))
            for (i in 2:length(colnames(extra_metadata))){  # skip the first column as it should be library_id
                current_column <- colnames(extra_metadata)[i]
                print(paste0("  Adding/replacing values in the '", current_column, "' column"))
                seurat_data[[current_column]] <- extra_metadata[[current_column]][match(seurat_data$new.ident, extra_metadata$library_id)]
            }
            return (seurat_data)
        } else {
            print(paste("Extra metadata loaded from ", args$condition, "is malformed. Skipping."))
            return (seurat_data)
        }
    } else {
        print("Extra metadata is not provided")
        return (seurat_data)
    }
}


get_combined_de_data <- function(seurat_data, design_formula, reduced_formula, args, assay="RNA", matrix_slot="counts"){
    print(paste0("Running DESeq2 using LRT for '", matrix_slot, "' slot from '", assay, "' assay."))
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    
    Idents(seurat_data) <- args$splitby
    pct_data <- FoldChange(
        seurat_data,
        ident.1=args$first,
        ident.2=args$second,
        slot=matrix_slot
    )
    Idents(seurat_data) <- "new.ident"
    alpha_min <- pmax(pct_data$pct.1, pct_data$pct.2)
    names(x=alpha_min) <- rownames(x=pct_data)
    features_high_pct <- names(x=which(x=alpha_min>=args$minpct))
    selected_features <- features_high_pct[!features_high_pct %in% args$exgenes]
    counts <- GetAssayData(seurat_data, assay=assay, slot=matrix_slot)    # set assay in case GetAssayData won't take the default value
    counts <- counts[selected_features, , drop=FALSE]
    pct_data <- pct_data[rownames(counts), , drop=FALSE]                  # to guarantee the order of rows is the same as in the counts
    col_data <- seurat_data@meta.data
    col_data[[args$timeby]] <- as.factor(col_data[[args$timeby]])         # double check if we really need it. DESeq complains if it's not factor
    col_data[[args$splitby]] <- as.factor(col_data[[args$splitby]])
    print("Time factor levels")
    print(levels(col_data[[args$timeby]]))
    print("Group factor levels")
    col_data[[args$splitby]] <- relevel(col_data[[args$splitby]], args$first)  # relevel so the direction of comparison will be args$first vs args$second
    print(levels(col_data[[args$splitby]]))
    deset_data <- DESeqDataSetFromMatrix(countData=counts, colData=col_data, design=design_formula)
    scran_size_factors <- scran::computeSumFactors(deset_data)            # setting size factor as recommended here https://github.com/mikelove/zinbwave-deseq2/blob/master/zinbwave-deseq2.knit.md
    sizeFactors(deset_data) <- sizeFactors(scran_size_factors)
    deset_data <- DESeq(
        deset_data,
        test="LRT",
        fitType="glmGamPoi",
        useT=TRUE,
        minmu=1e-6,
        minReplicatesForReplace=Inf,
        reduced=reduced_formula,
        quiet=FALSE,
        parallel=TRUE
    )
    DefaultAssay(seurat_data) <- backup_assay
    combined_data <- list(
        deset_data=deset_data,
        pct_data=pct_data[, c("pct.1", "pct.2"), drop=FALSE]
    )
    return (combined_data)
}


get_filtered_data <- function(seurat_data, args){
    if (!is.null(args$groupby) && !is.null(args$select)){
        print(
            paste0(
                "Filtering Seurat data to include only '",
                paste(args$select, collapse=" "), "' values from '",
                args$groupby, "' metadata column"
            )
        )
        Idents(seurat_data) <- args$groupby
        seurat_data <- subset(seurat_data, idents=args$select)
        Idents(seurat_data) <- "new.ident"
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
    DefaultAssay(seurat_data) <- assay
    Idents(seurat_data) <- args$splitby
    all_features <- as.vector(as.character(rownames(seurat_data)))
    selected_features <- all_features[!all_features %in% args$exgenes]
    diff_expr_markers <- FindMarkers(
        seurat_data,
        slot=slot,
        ident.1=args$first,
        ident.2=args$second,
        features=selected_features,
        logfc.threshold=args$minlogfc,
        min.pct=args$minpct,
        test.use=args$testuse,
        min.diff.pct=min_diff_pct,
        base=2,                      # to make sure we use log2 scale
        verbose=FALSE
    ) %>% rownames_to_column(var="gene")
    Idents(seurat_data) <- "new.ident"
    DefaultAssay(seurat_data) <- backup_assay
    return (diff_expr_markers)
}


export_counts_plot <- function(data, features, rootname, splitby, timeby, x_label, y_label, legend_title, plot_title, combine_guides=NULL, palette="Paired", alpha=0.1, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plots = list()
            for (i in 1:length(features)){
                current_feature <- features[i]
                normalized_counts <- plotCounts(
                    deseq_data,
                    current_feature, 
                    intgroup=c(timeby, splitby),
                    returnData=TRUE
                )
                normalized_counts[[timeby]] <- as.numeric(as.character(normalized_counts[[timeby]]))
                plots[[i]] <- ggplot(
                            normalized_counts,
                            aes_string(x=timeby, y="count", color=splitby, group=splitby)
                        ) + 
                        geom_jitter(alpha=alpha) + stat_summary(fun.y=mean, geom="line") +
                        ggtitle(current_feature) +
                        theme_gray() +
                        xlab(x_label) +
                        ylab(y_label) +
                        scale_y_log10()
            }
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export counts plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export counts plot to ", rootname, ".(png/pdf)", sep=""))
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
    parser <- ArgumentParser(description='Runs time course differential expression analysis for a subset of cells between two selected conditions')
    parser$add_argument(
        "--rds",
        help=paste(
            "Path to the RDS file to load Seurat object from.",
            "RDS file with the proper structure can be produced by run_seurat.R script."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--condition",
        help=paste(
            "Path to the TSV/CSV file to optionally extend Seurat object metadata. First",
            "column 'library_id' should include all unique values from the 'new.ident'",
            "column of the Seurat object metadata. All other columns will be added to the",
            "Seurat object metadata. If any of the provided in this file columns were",
            "already present in the Seurat object metadata, they will be overwritten.",
            "Default: no metadata columns will be added or overwritten"
        ),
        type="character"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Field from the Seurat object metadata to split cells into groups",
            "for --first vs --second differential expression analysis.",
            "May include values from the metadata fields added with --condition."
        ),
        type="character", required="True"
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
        "--timeby",
        help=paste(
            "Field from the Seurat object metadata to define time points.",
            "Should bo convertable to numeric values.",
            "May include values from the metadata fields added with --condition."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Field from the Seurat object metadata to group cells for optional",
            "subsetting (for example, clustering or predicted cell type field).",
            "May include values from the metadata fields added with --condition."
        ),
        type="character"
    )
    parser$add_argument(
        "--select",
        help=paste(
            "Value(s) from the column set with --groupby to optionally subset cells",
            "before running differential expression analysis.",
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
        "--exgenes",
        help=paste(
            "Genes to be excluded from the differential expression analysis.",
            "Default: include all genes"
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

print("Trying to extend Seurat object with extra metadata")
seurat_data <- extend_metadata(seurat_data, args)

design_formula <- as.formula(
    paste0("~", args$splitby, "+", args$timeby, "+", args$splitby, ":", args$timeby)
)
reduced_formula <- as.formula(
    paste0("~", args$splitby, "+", args$timeby)
)
contrast <- c(args$splitby, args$first, args$second)

print(
    paste0(
        "Identifying differentially expressed genes for time course experiment over '", args$timeby, "' column ",
        "between '", args$first, "' and '", args$second, "' conditions from '", args$groupby, "' column. ",
        ifelse(
            (!is.null(args$groupby) && !is.null(args$select)),
            paste0("Using cells only from '", paste(args$select, collapse=" "), "' groups based on '", args$groupby, "' metadata column. "),
            ""
        ),
        "Filtering result by (Padj <= ", args$maxpvadj, "). ",
        ifelse(
            is.null(args$genes),
            paste0("Highlighting ", 2*args$topn, " genes with the highest abs(log2FC)."),
            ""
        )
    )
)

seurat_data <- get_filtered_data(seurat_data, args)



combined_de_data <- get_combined_de_data(seurat_data, design_formula, reduced_formula, args)
deseq_data <- combined_de_data$deset_data
print(resultsNames(deseq_data))

deseq_results <- results(
    deseq_data,
    # contrast=contrast,
    contrast=list(c("condition_PP_vs_PV", "conditionPP.day4")),
    alpha=args$maxpvadj  # recommended to set to our FDR threshold https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
)
print(head(deseq_results))


diff_expr_genes <- deseq_results %>% data.frame() %>%
                   rownames_to_column(var="gene") %>%
                   filter(.$padj<=args$maxpvadj) %>%
                   arrange(desc(log2FoldChange))
diff_expr_genes <- cbind(
    diff_expr_genes,
    combined_de_data$pct_data[diff_expr_genes$gene, , drop=FALSE]
)

max_log2FC <- max(diff_expr_genes$log2FoldChange[is.finite(diff_expr_genes$log2FoldChange)])
min_log2FC <- min(diff_expr_genes$log2FoldChange[is.finite(diff_expr_genes$log2FoldChange)])
diff_expr_genes <- diff_expr_genes %>%                                                        # Maybe we don't need ot. Replaces all Inf and -Inf with numbers
                   mutate("log2FoldChange"=ifelse(
                       .$log2FoldChange==Inf, max_log2FC+0.05*abs(max_log2FC), ifelse(
                           .$log2FoldChange==-Inf, min_log2FC-0.05*abs(min_log2FC), .$log2FoldChange)
                       )
                   )

topn_diff_expr_genes <- diff_expr_genes %>% filter(row_number() > max(row_number()) - all_of(args$topn) | row_number() <= all_of(args$topn))
highlight_genes <- as.vector(as.character(topn_diff_expr_genes[, "gene"]))  # default genes to highlight
if (!is.null(args$genes)){
    print("Check genes of interest to include only those that are differentially expressed")
    args$genes <- unique(args$genes)
    args$genes <- args$genes[args$genes %in% as.vector(as.character(diff_expr_genes[, "gene"]))]
    highlight_genes <- args$genes
}
print(paste("Genes to highlight", paste(highlight_genes, collapse=", ")))


export_counts_plot(
    data=deseq_data,
    features=highlight_genes,
    rootname=paste(args$output, "norm_counts", sep="_"),
    splitby=args$split,
    timeby=args$timeby,
    x_label="Time",
    y_label="Counts",
    legend_title=deparse(substitute(args$split)),
    plot_title="Normalized counts",
    combine_guides="collect",
    width=5000, height=5000, resolution=200,
    pdf=args$pdf
)

export_volcano_plot(
    data=diff_expr_genes,
    rootname=paste(args$output, "diff_expr_genes", sep="_"),
    x_axis="log2FoldChange",
    y_axis="padj",
    x_cutoff=0,
    y_cutoff=args$maxpvadj,
    x_label=bquote(~Log[2] ~ "FC"),
    y_label=bquote(~-Log[10] ~ italic(Padj)),
    plot_title="Differentially expressed genes",
    plot_subtitle=paste0(
        args$first, " vs ", args$second, " for cells split by ", args$splitby,
        " (Padj <= ", args$maxpvadj, ")",
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






quit(save="no", status=0, runLast=FALSE)
















# Filter DESeq2 output
res_filtered <- as.data.frame(deseq_res[,c(2,5,6)])
res_filtered$log2FoldChange[is.na(res_filtered$log2FoldChange)] = 0;
res_filtered[is.na(res_filtered)] = 1;

# Export results to TSV file
expression_data_df = res_filtered
expression_data_df[,"-LOG10(pval)"] <- -log(as.numeric(expression_data_df$pval), 10)
expression_data_df[,"-LOG10(padj)"] <- -log(as.numeric(expression_data_df$padj), 10)





 
