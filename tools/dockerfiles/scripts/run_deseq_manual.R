#!/usr/bin/env Rscript
options(warn=-1)
options("width"=400)

suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(argparse))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(BiocParallel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(EnhancedVolcano))


# Useful links
# 1. https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md


COUNTS_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")


set_threads <- function (threads) {
    register(MulticoreParam(threads))
}


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


load_metadata <- function(args){
    metadata <- read.table(
        args$metadata,
        sep=get_file_type(args$metadata),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )  %>% remove_rownames() %>% column_to_rownames("sample") %>% mutate_at(colnames(.), factor)
    if (!is.null(args$base)){
        print(
            paste(
                "Attempting to relevel metadata table based on",
                paste(args$base, collapse=", "), "values"
            )
        )
        for (i in 1:length(args$base)) {
            current_base_level <- args$base[i]
            current_column <- colnames(metadata)[i]
            tryCatch(
                expr = {
                    metadata[[current_column]] <- relevel(metadata[[current_column]], current_base_level)
                    print(paste("Setting", current_base_level, "as a base level for", current_column))
                },
                error = function(e){
                    print(paste("Failed to set", current_base_level, "as a base level for", current_column))
                }
            )
        }
    }
    return (metadata)
}


load_expression_data <- function(args, counts_colname=COUNTS_COL, rpkm_colname=RPKM_COL, intersect_by=INTERSECT_BY) {
    collected_expression_data <- NULL
    for (i in 1:length(args$expression)) {
        location <- args$expression[i]
        alias <- args$aliases[i]
        expression_data <- read.table(location, sep=get_file_type(location), header=TRUE, stringsAsFactors=FALSE)
        print(paste("Loading", nrow(expression_data), "rows from", location, "as", alias))
        colnames(expression_data)[colnames(expression_data) == counts_colname] <- paste(alias, counts_colname, sep=" ")
        colnames(expression_data)[colnames(expression_data) == rpkm_colname] <- paste(alias, rpkm_colname, sep=" ")
        if (is.null(collected_expression_data)){
            collected_expression_data <- expression_data
        } else {
            collected_expression_data <- merge(collected_expression_data, expression_data, by=intersect_by, sort = FALSE)
        }
    }
    print(paste("Number of rows common for all loaded files ", nrow(collected_expression_data), sep=""))
    collected_expression_data <- collected_expression_data %>% remove_rownames()
    if (args$ftype == "gene"){
        collected_expression_data <- collected_expression_data %>% column_to_rownames("GeneId")
    } else {
        collected_expression_data <- collected_expression_data %>% column_to_rownames("RefseqId")
    }
    all_features <- as.vector(as.character(rownames(collected_expression_data)))
    selected_features <- all_features[!all_features %in% args$exfeatures]
    excluded_features <- all_features[all_features %in% args$exfeatures]   # not used elsewhere, only to print to the console
    if (length(excluded_features) > 0){
        print(paste("Following features will be excluded from the differential expression analysis:", paste(excluded_features, collapse=", ")))
    }
    collected_expression_data <- collected_expression_data[selected_features,]

    print(paste("Keeping only those features which total counts for all samples are bigger then", args$mincounts))
    counts_columns <- grep(
        paste(COUNTS_COL, sep=""),
        colnames(collected_expression_data),
        value=TRUE,
        ignore.case=TRUE
    )
    keep <- rowSums(collected_expression_data[counts_columns]) > args$mincounts
    collected_expression_data <- collected_expression_data[keep,]
    print(paste("Number of the remaining features: ", nrow(collected_expression_data), sep=""))
    return (collected_expression_data)
}


assert_args <- function(args){
    if (length(args$expression) != length(args$aliases)){
            print("Exiting: --expression and --aliases have different number of values")
            quit(save = "no", status = 1, runLast = FALSE)
    }
    return (args)
}


get_highlight_features <- function(diff_expr_features, args){
    if (args$usepvalue){
        print(paste("Filtering DESeq results to include only features with pvalue <= ", args$padj))
        filt_diff_expr_features <- diff_expr_features %>% filter(.$pvalue<=args$padj)
    } else {
        print(paste("Filtering DESeq results to include only features with padj <= ", args$padj))
        filt_diff_expr_features <- diff_expr_features %>% filter(.$padj<=args$padj)
    }
    filt_diff_expr_features <- filt_diff_expr_features %>% arrange(desc(log2FoldChange))

    print(paste("Number of significantly differentially expressed features:", nrow(filt_diff_expr_features)))

    topn_diff_expr_features <- filt_diff_expr_features %>% filter(row_number() > max(row_number()) - all_of(args$topn) | row_number() <= all_of(args$topn))
    highlight_features <- as.vector(as.character(topn_diff_expr_features[, "feature"]))                      # default features to highlight
    if (!is.null(args$features)){
        print("Check features of interest to include only those that are differentially expressed regardless of significance tresholds")
        args$features <- unique(args$features)
        args$features <- args$features[args$features %in% as.vector(as.character(diff_expr_features[, "feature"]))]
        highlight_features <- args$features
    }
    print(paste("Features to highlight", paste(highlight_features, collapse=", ")))
    return (highlight_features)
}


get_contrast <- function(design_formula, deseq_data, metadata, args){
    model <- model.matrix(design_formula, metadata)
    print("Model matrix")
    print(model)
    for (i in 1:length(colnames(metadata))){
        current_column <- colnames(metadata)[i]
        unique_keys <- unique(metadata[[current_column]])
        for (j in 1:length(unique_keys)){
            key <- as.character(unique_keys[j])
            print(paste("Subsetting model matrix for", current_column, "column", "with value", key))
            subset <- colMeans(model[deseq_data[[current_column]] == key, ])
            print(subset)
            assign(key, subset)
        }
    }
    print(paste("Evaluating contrast", args$contrast))
    contrast <- eval(parse(text=args$contrast))
    print(contrast)
    return (contrast)
}


get_diff_expr_data <- function(expression_data, metadata, args){
    print("Loading design formula")
    design_formula <- as.formula(args$design)
    if (!args$wald){
        print("Loading reduced formula")
        reduced_formula <- as.formula(args$reduced)
    }

    print("Selecting all columns with read counts data")
    counts_columns <- grep(
        paste(COUNTS_COL, sep=""),
        colnames(expression_data),
        value=TRUE,
        ignore.case=TRUE
    )
    counts_data = expression_data[counts_columns]
    colnames(counts_data) <- lapply(
        colnames(counts_data),
        function(s){
            paste(head(unlist(strsplit(s, " ", fixed=TRUE)), -1), collapse=" ")
        }
    )
    print("Reordering read counts data columns based on the row names from the provided metadata file")
    counts_data <- counts_data[, rownames(metadata)]

    if ( all(colnames(counts_data) != rownames(metadata)) ){   # safety measure
        print("Assert failed: columns order of the counts data should be the same as rows order in metadata")
        quit(save = "no", status = 1, runLast = FALSE)
    }

    print("Loadind data to DESeq2")
    deseq_data <- DESeqDataSetFromMatrix(
        countData=counts_data,
        colData=metadata,
        design=design_formula
    )

    if (args$wald){
        print("Forced to use Wald test to calculate p-values instead of LRT")
        deseq_data <- DESeq(
            deseq_data,
            quiet=TRUE,
            parallel=TRUE,
            BPPARAM=MulticoreParam(args$threads)  # add it here as well just in case
        )
    } else {
        print("Using LRT test to calculate p-values")
        deseq_data <- DESeq(
            deseq_data,
            test="LRT",
            reduced=reduced_formula,
            quiet=TRUE,
            parallel=TRUE,
            BPPARAM=MulticoreParam(args$threads)  # add it here as well just in case
        )
    }
    print("Estimated effects")
    print(resultsNames(deseq_data))
    deseq_results <- results(
        deseq_data,
        contrast=get_contrast(design_formula, deseq_data, metadata, args),
        parallel=TRUE,
        BPPARAM=MulticoreParam(args$threads)  # add it here as well just in case
    )

    print("Results description")
    print(mcols(deseq_results))
    print(head(deseq_results))

    print("Adding extra columns to the DESeq output")
    rpkm_columns <- grep(              # for column ordering only
        paste(RPKM_COL, sep=""),
        colnames(expression_data),
        value=TRUE,
        ignore.case=TRUE
    )
    deseq_results <- expression_data %>%
                     bind_cols(as.data.frame(deseq_results)) %>%
                     na.omit() %>%                                           # exclude all rows where NA is found in any column (comes from DESeq)
                     rownames_to_column(var="feature") %>%                   # will make "feature" the first columns instead of GeneId or RefseqId depending on --ftype value
                     relocate(any_of(rpkm_columns), .after=last_col()) %>%   # move all rpkm columns to the end
                     relocate(any_of(counts_columns), .after=last_col())     # move all read counts columns to the end
    print(paste("Number of differentially expressed features after excluding NA:", nrow(deseq_results)))
    return (list(features=deseq_results, raw=deseq_data))
}


export_volcano_plot <- function(data, rootname, x_axis, y_axis, x_cutoff, y_cutoff, x_label, y_label, plot_title, plot_subtitle, caption, features=NULL, label_column="feature", pdf=FALSE, width=1200, height=800, resolution=100){
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
                        selectLab=features,
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
                    theme(legend.position="none", plot.subtitle=element_text(size=8, face="italic", color="gray30"))

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


export_counts_plot <- function(data, metadata, features, captions, rootname, splitby, groupby, x_label, y_label, legend_title, plot_title, plot_subtitle, combine_guides=NULL, alpha=1, pdf=FALSE, width=1600, height=1400, resolution=100){
    tryCatch(
        expr = {
            plots = list()
            normalized_counts <- t(as.data.frame(assay(rlog(data, blind=FALSE))))    # transposed rlog normalized counts
            normalized_counts <- normalized_counts[, features, drop=FALSE]           # to include only features we need
            normalized_counts <- merge(normalized_counts, metadata, by=0) %>%        # merge by row names
                                 remove_rownames() %>%
                                 column_to_rownames("Row.names")
            for (i in 1:length(features)){
                current_feature <- features[i]
                current_caption <- captions[i]
                counts <- normalized_counts[, c(current_feature, colnames(metadata)), drop=FALSE] %>%
                          dplyr::rename("count"=all_of(current_feature))
                plots[[i]] <- ggplot(
                                  counts,
                                  aes_string(x=groupby, y="count", color=splitby, group=splitby)
                              ) +
                              geom_point(alpha=alpha) +
                              stat_summary(fun=mean, geom="line") +
                              ggtitle(current_feature, subtitle=current_caption) +
                              theme_gray() +
                              xlab(x_label) +
                              ylab(y_label) +
                              theme(plot.subtitle=element_text(size=8, face="italic", color="gray30"))
            }
            combined_plots <- wrap_plots(plots, guides=combine_guides) +
                              plot_annotation(
                                  title=plot_title,
                                  subtitle=plot_subtitle,
                                  theme=theme(plot.subtitle=element_text(size=8, face="italic", color="gray30"))
                              )

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


export_pca_plot <- function(data, rootname, intgroup, plot_title, plot_subtitle, ntop=500, palette="Paired", alpha=0.5, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            rlog_data <- rlog(data, blind=TRUE)     # rlog normalized counts
            pca_data <- plotPCA(
                rlog_data,
                intgroup=intgroup,
                ntop=ntop,               # ntop the most variable genes
                returnData=TRUE
            )
            percentVar <- round(100 * attr(pca_data, "percentVar"))
            plot <- ggplot(pca_data, aes(PC1, PC2, color=group)) +
                    geom_point(size=4, shape=19, alpha=alpha) +
                    xlab(paste0("PC1: ",percentVar[1], "% variance")) +
                    ylab(paste0("PC2: ",percentVar[2], "% variance")) + 
                    ggtitle(plot_title, subtitle=plot_subtitle) +
                    geom_text_repel(
                        aes(label=name),
                        point.padding=0.5,
                        box.padding=0.5,
                        check_overlap=TRUE,
                        show.legend=FALSE
                    ) +
                    theme(plot.subtitle=element_text(size=8, face="italic", color="gray30"))
                    scale_fill_brewer(palette=palette)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export PCA plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export PCA plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_data <- function(data, location, row_names=FALSE, col_names=TRUE, quote=FALSE, digits=NULL){
    tryCatch(
        expr = {
            if (!is.null(digits)){
                data <- format(data, digits=digits)
            }
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
    parser <- ArgumentParser(description="Run DeSeq2 with manual control over major parameters")
    parser$add_argument(
        "--expression",
        help=paste(
            "Path to the TSV/CSV files with expression data. All files should have the following header: ",
            "RefseqId GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm"
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--aliases",
        help=paste(
            "Unique names for files provided in --expression, no special characters or spaces are allowed.",
            "Number and order of the names should corresponds to values from --expression"
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to provide metadata for the samples from --expression.",
            "First column should have the name 'sample', other columns may have arbitrary names.",
            "The values from the 'sample' column should correspond to the values provided in --aliases.",
            "For a proper --contrast intepretation, values defined in each column should not be used in others."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--design",
        help=paste(
            "Design formula. Should start with ~ and include terms from the --metadata table"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--reduced",
        help=paste(
            "Reduced formula to compare against with the term(s) of interest removed.",
            "Should start with ~. Ignored when run with --wald"
        ),
        type="character"
    )
    parser$add_argument(                                                                         # https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
        "--contrast",                                                                            # https://github.com/tavareshugo/tutorial_DESeq2_contrasts
        help=paste(
            "Contrast to be be applied for the output, formatted as a mathematical formula",
            "of values from the --metadata table"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--base",
        help=paste(
            "Value from each column of metadata file to be set as base levels.",
            "Number and order of provided values should correspond the order of columns",
            "in --metadata file.",
            "Default: define base levels alphabetically for each columns of metadata table"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--ftype",
        help=paste(
            "Feature type to use for differential expression.",
            "If set to 'gene', use 'GeneId' column from the provided in --expression files.",
            "If set to 'transcript', use 'RefseqId' from the provided in --expression files.",
            "Default: gene"
        ),
        type="character", default="gene",
        choices=c("gene", "transcript")
    )
    parser$add_argument(
        "--mincounts",
        help=paste(
            "Keep only those features where the total number of counts for all samples",
            "is bigger than this value.",
            "Default: 0"
        ),
        type="integer", default=0
    )
    parser$add_argument(
        "--wald",
        help=paste(
            "Use pair-wise Wald test instead of LRT. --reduced parameter will be ignored",
            "Default: use LRT test"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Used only in plots. Column from the metadata file to split samples into categories.",
            "Default: the first after the 'sample' column from the metadata file"
        ),
        type="character"
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Used only in plots. Column from the metadata file to combine samples into groups.",
            "Default: the last column from the metadata file"
        ),
        type="character"
    )
    parser$add_argument(
        "--features",
        help=paste(
            "Used only in plots. Features of interest to label on the generated plots.",
            "Default: --topn N features with the highest and the lowest log2 fold change",
            "expression values."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--exfeatures",
        help=paste(
            "Used only in plots. Features to be excluded from the differential expression analysis.",
            "Default: include all features"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--topn",
        help=paste(
            "Used only in plots. Show N features with the highest and N features with the lowest log2 fold",
            "change expression values. Ignored with --features.",
            "Default: 10"
        ),
        type="integer", default=10
    )
    parser$add_argument(
        "--padj",
        help=paste(
            "Used only in plots. Output only features with adjusted P-value not bigger than this treshold.",
            "Default: 0.05"
        ),
        type="double", default=0.05
    )
    parser$add_argument(
        "--usepvalue",
        help=paste(
            "Used only in plots. Treat --padj as a theshold for P-value",
            "Default: --padj defines the treshold for adjusted P-value"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help=paste(
            "Output prefix for generated files"
        ),
        type="character", default="./deseq"
    )
    parser$add_argument(
        "--threads",
        help=paste(
            "Threads number"
        ),
        type="integer", default=1
    )
    args <- assert_args(parser$parse_args(commandArgs(trailingOnly = TRUE)))
    return (args)
}


# Parse arguments
args <- get_args()

print("Used parameters")
print(args)

print(paste("Setting parallelizations threads to", args$threads))
set_threads(args$threads)

print("Loading expression data")
expression_data <- load_expression_data(args)

print(paste("Loading metadata from", args$metadata))
metadata <- load_metadata(args)
print(metadata)
print(str(metadata))

print("Running differential expression analysis")
diff_expr_data <- get_diff_expr_data(expression_data, metadata, args)

print("Identifying features to highlight")
highlight_features <- get_highlight_features(diff_expr_data$features, args)

export_volcano_plot(
    data=diff_expr_data$features,                                   # this is not filtered differentially expressed features
    rootname=paste(args$output, "volcano_plot", sep="_"),
    x_axis="log2FoldChange",
    y_axis=ifelse(args$usepvalue, "pvalue", "padj"),
    x_cutoff=0,
    y_cutoff=args$padj,
    x_label="log2FoldChange",
    y_label=ifelse(args$usepvalue, "-log10 Pvalue", "-log10 Padj"),
    plot_title="Differentially expressed features",
    plot_subtitle=paste0(
        "Differentially expressed features with ",
        ifelse(args$usepvalue, "pvalue", "padj"), "<=", args$padj,
        " for contrast ", args$contrast
    ),
    caption=paste(nrow(diff_expr_data$features), "features"),
    features=highlight_features,
    pdf=args$pdf
)

captions <- diff_expr_data$features %>%                                           # this is not filtered differentially expressed features
            filter(.$feature %in% highlight_features) %>%
            mutate(
                "caption" = paste0(
                    "log2FC=", round(.$log2FoldChange, 5), "\npvalue=", round(.$pvalue, 5), "\npadj=", round(.$padj, 5)
                )
            )
captions <- captions[match(highlight_features, captions$feature), ] %>%           # to guarantee that the order of captions will correspond to the order of highlight_genes
            pull(caption)

args$splitby <- ifelse(is.null(args$splitby), colnames(metadata)[1], args$splitby)
args$groupby <- ifelse(is.null(args$groupby), colnames(metadata)[-1], args$groupby)

export_counts_plot(
    data=diff_expr_data$raw,
    metadata=metadata,
    features=highlight_features,
    captions=captions,
    rootname=paste(args$output, "counts_plot", sep="_"),
    splitby=args$splitby,
    groupby=args$groupby,
    x_label=args$groupby,
    y_label="Counts",
    legend_title=deparse(substitute(args$splitby)),
    plot_title="rlog-normalized counts",
    plot_subtitle=paste(
        "Split by", args$splitby,
        "grouped by", args$groupby,
        "for contrast", args$contrast
    ),
    combine_guides="collect",
    pdf=args$pdf
)

export_pca_plot(
    data=diff_expr_data$raw,
    rootname=paste(args$output, "pca_plot", sep="_"),
    intgroup=colnames(metadata),
    plot_title="PCA plot of rlog-normalized counts",
    plot_subtitle=paste(
        "Based on the top", 500, "genes selected by the highest row variance"
    ),
    ntop=500,
    pdf=args$pdf
)

print("Exporting differentially expressed genes")
export_data(
    diff_expr_data$features,                                        # this is not filtered differentially expressed features
    location=paste(args$output, "diff_expr_genes.tsv", sep="_"),
    digits=5
)
