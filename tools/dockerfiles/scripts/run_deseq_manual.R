#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(cmapR))
suppressMessages(library(dplyr))
suppressMessages(library(limma))
suppressMessages(library(DESeq2))
suppressMessages(library(hopach))
suppressMessages(library(Glimma))
suppressMessages(library(argparse))
suppressMessages(library(tidyverse))
suppressMessages(library(patchwork))
suppressMessages(library(htmlwidgets))
suppressMessages(library(BiocParallel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(SummarizedExperiment))


# Useful links
# 1. https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
# 2. https://www.statlect.com/glossary/design-matrix#:~:text=A%20design%20matrix%20is%20a,each%20column%20to%20a%20characteristic.
# 3. https://paasp.net/accurate-design-of-in-vitro-experiments-why-does-it-matter/
# 4. https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
# 5. https://github.com/tavareshugo/tutorial_DESeq2_contrasts


COUNTS_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")


set_cpus <- function (cpus) {
    register(MulticoreParam(cpus))
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
            tryCatch(expr={dev.off()}, error=function(err){print(paste("Called  dev.off() with error -", err))})
            print(paste0("Failed to export volcano plot to ", rootname, ".(png/pdf) with error - ", e))
        }
    )
}


export_pca_plot <- function(norm_counts_data, rootname, intgroup, plot_title, plot_subtitle, ntop=500, palette="Paired", alpha=0.5, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            pca_data <- plotPCA(
                norm_counts_data,
                intgroup=intgroup,
                ntop=ntop,               # ntop the most variable features
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
            tryCatch(expr={dev.off()}, error=function(err){print(paste("Called  dev.off() with error -", err))})
            print(paste0("Failed to export PCA plot to ", rootname, ".(png/pdf) with error - ", e))
        }
    )
}


export_mds_html_plot <- function(norm_counts_data, location){
    tryCatch(
        expr = {
            htmlwidgets::saveWidget(
                glimmaMDS(
                    x=assay(norm_counts_data),
                    groups=as.data.frame(SummarizedExperiment::colData(norm_counts_data)),
                    labels=rownames(SummarizedExperiment::colData(norm_counts_data))
                ),
                file=location
            )
        },
        error = function(e){
            print(paste0("Failed to export MDS plot to ", location, " with error - ", e))
        }
    )
}


get_highlight_features <- function(diff_expr_features, args){
    print(paste("Filtering DESeq results to include only features with padj <= ", args$padj))
    filt_diff_expr_features <- diff_expr_features %>% filter(.$padj<=args$padj)
    filt_diff_expr_features <- filt_diff_expr_features %>% arrange(desc(log2FoldChange))

    print(paste("Number of significantly differentially expressed features:", nrow(filt_diff_expr_features)))

    topn_diff_expr_features <- filt_diff_expr_features %>% filter(row_number() > max(row_number()) - 10 | row_number() <= 10)
    highlight_features <- as.vector(as.character(topn_diff_expr_features[, "feature"]))                      # default features to highlight
    if (!is.null(args$label)){
        print("Check features of interest to include only those that are differentially expressed regardless of significance tresholds")
        args$label <- unique(args$label)
        args$label <- args$label[args$label %in% as.vector(as.character(diff_expr_features[, "feature"]))]
        highlight_features <- args$label
    }
    print(paste("Features to highlight", paste(highlight_features, collapse=", ")))
    return (highlight_features)
}


export_plots <- function(diff_expr_data, norm_counts_data, metadata, args){
    print("Identifying features to highlight")
    highlight_features <- get_highlight_features(diff_expr_data$res, args)
    export_volcano_plot(
        data=diff_expr_data$res,                                   # this is not filtered differentially expressed features
        rootname=paste(args$output, "volcano_plot", sep="_"),
        x_axis="log2FoldChange",
        y_axis="padj",
        x_cutoff=0,
        y_cutoff=args$padj,
        x_label="log2FoldChange",
        y_label="-log10 Padj",
        plot_title="Differentially expressed features",
        plot_subtitle=paste0(
            "Differentially expressed features with padj", "<=", args$padj,
            ifelse(args$contrast, paste0(" for contrast ", args$contrast), "")
        ),
        caption=paste(nrow(diff_expr_data$res), "features"),
        features=highlight_features,
        pdf=args$pdf
    )
    export_pca_plot(
        norm_counts_data=norm_counts_data,
        rootname=paste(args$output, "pca_plot", sep="_"),
        intgroup=colnames(metadata),
        plot_title=paste0("PCA plot of ", args$norm, "-normalized read counts"),
        plot_subtitle=paste(
            "Based on the top 500 features.",
            ifelse(!is.null(args$remove), paste("Remove the effect of", args$remove), "")
        ),
        ntop=500,
        pdf=args$pdf
    )
    export_mds_html_plot(
        norm_counts_data=norm_counts_data,
        location=paste(args$output, "mds_plot.html", sep="_")
    )
}


export_gct <- function(counts_mat, row_metadata, col_metadata, location){
    tryCatch(
        expr = {
            row_metadata <- row_metadata %>% rownames_to_column("id") %>% mutate_at("id", as.vector)
            col_metadata <- col_metadata %>% rownames_to_column("id") %>% mutate_at("id", as.vector)
            gct_data <- new(
                "GCT",
                mat=counts_mat[row_metadata$id, col_metadata$id],       # to guarantee the order and number of row/columns
                rdesc=row_metadata,
                cdesc=col_metadata
            )
            write_gct(
                ds=gct_data,
                ofile=location,
                appenddim=FALSE
            )
            print(paste("Exporting GCT data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export GCT data to", location, sep=" "))
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
            tryCatch(
                expr = {
                    current_column <- colnames(metadata)[i]
                    metadata[[current_column]] <- relevel(metadata[[current_column]], current_base_level)
                    print(paste("Setting", current_base_level, "as a base level for", current_column))
                },
                error = function(e){
                    print(paste("Failed to set", current_base_level, "as a base level"))
                }
            )
        }
    }
    print("Loaded metadata")
    print(metadata)
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
    if (args$type == "gene"){
        collected_expression_data <- collected_expression_data %>% column_to_rownames("GeneId")
    } else {
        collected_expression_data <- collected_expression_data %>% column_to_rownames("RefseqId")
    }
    all_features <- as.vector(as.character(rownames(collected_expression_data)))
    selected_features <- all_features[!all_features %in% args$exclude]
    excluded_features <- all_features[all_features %in% args$exclude]   # not used elsewhere, only to print to the console
    if (length(excluded_features) > 0){
        print(paste("Following features will be excluded from the differential expression analysis:", paste(excluded_features, collapse=", ")))
    }
    collected_expression_data <- collected_expression_data[selected_features, ]

    print("Keeping only those features which total counts for all samples are bigger then 0")
    counts_columns <- grep(
        paste(COUNTS_COL, sep=""),
        colnames(collected_expression_data),
        value=TRUE,
        ignore.case=TRUE
    )
    keep <- rowSums(collected_expression_data[counts_columns]) > 0
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

    print("Selecting all columns with raw read counts data.")
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
    print(
        paste(
            "Reordering read counts data columns based on the",
            "row names from the provided metadata file."
        )
    )
    counts_data <- counts_data[, as.vector(rownames(metadata))]

    if ( all(colnames(counts_data) != rownames(metadata)) ){              # safety measure
        print(
            paste(
                "Assert failed: columns order of the counts data",
                "should be the same as rows order in metadata."
            )
        )
        quit(save = "no", status = 1, runLast = FALSE)
    }

    print("Loadind data to DESeq2")
    deseq_data <- DESeqDataSetFromMatrix(
        countData=counts_data,
        colData=metadata,
        design=as.formula(args$design)
    )

    if (!is.null(args$reduced)){
        print("Using LRT test to calculate p-values")
        deseq_data <- DESeq(
            deseq_data,
            test="LRT",
            reduced=as.formula(args$reduced),
            quiet=TRUE,
            parallel=TRUE,
            BPPARAM=MulticoreParam(args$cpus)  # add it here as well just in case
        )
    } else {
        print("Using Wald test to calculate p-values")
        deseq_data <- DESeq(
            deseq_data,
            quiet=TRUE,
            parallel=TRUE,
            BPPARAM=MulticoreParam(args$cpus)  # add it here as well just in case
        )
    }
    print("Estimated effects")
    print(resultsNames(deseq_data))

    if (is.null(args$contrast)){                         # can't fake missing contrast parameter so need this if statement
        print(paste(
            "Contrast is not provided. The last term",
            "from the design formula will be used."
        ))
        deseq_results <- results(
            deseq_data,
            parallel=TRUE,
            BPPARAM=MulticoreParam(args$cpus)         # add it here as well just in case
        )
    } else {
        print("Using user-provided contrast.")
        deseq_results <- results(
            deseq_data,
            contrast=get_contrast(as.formula(args$design), deseq_data, metadata, args),
            parallel=TRUE,
            BPPARAM=MulticoreParam(args$cpus)         # add it here as well just in case
        )
    }

    print("Results description")
    print(mcols(deseq_results))
    print(head(deseq_results))

    print("Adding extra columns to the DESeq output")
    rpkm_columns <- grep(                                                    # for column ordering only
        paste(RPKM_COL, sep=""),
        colnames(expression_data),
        value=TRUE,
        ignore.case=TRUE
    )
    deseq_results <- expression_data %>%
                     bind_cols(as.data.frame(deseq_results)) %>%
                     na.omit() %>%                                           # exclude all rows where NA is found in any column (comes from DESeq)
                     rownames_to_column(var="feature") %>%                   # will make "feature" the first columns instead of GeneId or RefseqId depending on --type value
                     relocate(any_of(rpkm_columns), .after=last_col()) %>%   # move all rpkm columns to the end
                     relocate(any_of(counts_columns), .after=last_col())     # move all read counts columns to the end
    print(paste("Number of differentially expressed features after excluding NA:", nrow(deseq_results)))
    print("DESeq2 results")
    print(head(deseq_results))
    return (list(
        res=deseq_results,
        raw=deseq_data
    ))
}


get_norm_counts_data <- function(deseq_data, metadata, args){
    if (args$norm == "vst"){
        print("Applying vst transformation (not blind to the experimental design)")
        norm_counts_data <- DESeq2::vst(deseq_data, blind=FALSE)
    } else {
        print("Applying rlog transformation (not blind to the experimental design)")
        norm_counts_data <- DESeq2::rlog(deseq_data, blind=FALSE)
    }
    if(!is.null(args$remove)){
        print(
            paste("Removing the effect of", args$remove, "from the normalized counts")
        )
        keep_formula <- paste(
            grep(
                args$remove,
                unlist(strsplit(args$design, "\\+")),
                value=TRUE,
                ignore.case=TRUE,
                invert=TRUE
            ),
            collapse="+"
        )
        print(paste("Formula to include conditions to be preserved", keep_formula))
        assay(norm_counts_data) <- limma::removeBatchEffect(
            assay(norm_counts_data),
            batch=norm_counts_data[[args$remove]],
            design=stats::model.matrix(stats::as.formula(keep_formula), metadata)
        )
    }
    print("Normalized read counts")
    print(head(assay(norm_counts_data)))
    print(dim(assay(norm_counts_data)))
    print(SummarizedExperiment::colData(norm_counts_data))
    return (norm_counts_data)
}


get_clustered_data <- function(expression_data, center, dist, transpose) {

    if (transpose){
        print("Transposing expression data")
        expression_data = t(expression_data)
    }
    if (!is.null(center)) {
        print(paste("Centering expression data by ", center, sep=""))
        if (center == "mean"){
            expression_data = expression_data - rowMeans(expression_data)    
        } else {
            expression_data = expression_data - rowMedians(data.matrix(expression_data))    
        }
    }
    print("Creating distance matrix")
    distance_matrix <- distancematrix(expression_data, dist)
    print("Running HOPACH")
    hopach_results <- hopach(expression_data, dmat=distance_matrix)

    print("Reordering expression data")
    expression_data <- expression_data[as.vector(hopach_results$final$order),]
    if (transpose){
        print("Transposing expression data")
        expression_data = t(expression_data)
    }

    return (list(order=as.vector(hopach_results$final$order), expression=expression_data))
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
    parser <- ArgumentParser(description="DESeq2 Multi-factor Analysis")
    parser$add_argument(
        "--expression",
        help=paste(
            "Path to the TSV/CSV files with expression data.",
            "All files should have the following header:",
            "RefseqId GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm"
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--aliases",
        help=paste(
            "Unique names for files provided in --expression,",
            "no special characters or spaces are allowed.",
            "Number and order of the names should corresponds",
            "to values from --expression."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to provide metadata for the",
            "samples from --expression. First column should have",
            "the name 'sample', other columns may have arbitrary names.",
            "The values from the 'sample' column should correspond to",
            "the values provided in --aliases. For a proper --contrast",
            "intepretation, values defined in each column should not be",
            "used in other columns. All metadata columns are treated as",
            "factors (no covariates are supported)."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--design",
        help=paste(
            "Design formula. Should start with ~ and include terms from",
            "the --metadata table."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--reduced",
        help=paste(
            "Reduced formula with the term(s) of interest removed.",
            "Should start with ~. If provided, force DESeq2 to run",
            "LRT test instead of the Wald."
        ),
        type="character"
    )
    parser$add_argument(
        "--contrast",
        help=paste(
            "Contrast to be be applied for the output, formatted as",
            "a mathematical formula of values from the --metadata table.",
            "If not provided, the last term from the design formula will",
            "be used."
        ),
        type="character"
    )
    parser$add_argument(
        "--base",
        help=paste(
            "Value(s) from each metadata file column(s) to be set as",
            "the base level(s). Number and order of provided values should",
            "correspond the order of columns in --metadata file. Default:",
            "define base levels alphabetically for each metadata column."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--type",
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
        "--exclude",
        help=paste(
            "Features to be excluded from the differential expression analysis.",
            "Default: include all features"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--norm",
        help=paste(
            "Read counts normalization for the exploratory visualization analysis.",
            "Use 'vst' for medium-to-large datasets (n > 30) and 'rlog' for",
            "small datasets (n < 30), when there is a wide range of sequencing",
            "depth across samples.",
            "Default: vst"
        ),
        type="character", default="vst",
        choices=c("vst", "rlog")
    )
    parser$add_argument(
        "--remove",
        help=paste(
            "Column from the metadata file to remove batch effect when",
            "exporting feature counts. All components that include this",
            "term will be removed from the design formula when correcting",
            "for batch effect. Default: do not remove batch effect from",
            "the exported counts"
        ),
        type="character"
    )
    parser$add_argument(
        "--cluster",
        help=paste(
            "Hopach clustering method to be run on normalized read counts for the",
            "exploratory visualization analysis. Default: do not run clustering"
        ),
        type="character",
        choices=c("row", "column", "both")
    )
    parser$add_argument(
        "--center",
        help=paste(
            "Apply mean centering for feature expression prior to running",
            "clustering by row. Ignored when --cluster is not row or both.",
            "Default: do not centered"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--label",
        help=paste(
            "Features of interest to label on the generated volcanot plot. Default:",
            "top 10 features with the highest and the lowest log2 fold change",
            "expression values."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--padj",
        help=paste(
            "In the exploratory visualization analysis output only features with",
            "adjusted P-value not bigger than this value. Default: 0.05"
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
        help=paste(
            "Output prefix for generated files"
        ),
        type="character", default="./deseq"
    )
    parser$add_argument(
        "--cpus",
        help="Number of cores/cpus to use. Default: 1",
        type="integer", default=1
    )
    args <- assert_args(parser$parse_args(commandArgs(trailingOnly = TRUE)))
    return (args)
}


# Parse arguments
args <- get_args()

print("Used parameters")
print(args)

print(paste("Setting parallelizations to", args$cpus, "cores"))
set_cpus(args$cpus)

print("Loading expression data")
expression_data <- load_expression_data(args)

print(paste("Loading metadata from", args$metadata))
metadata <- load_metadata(args)

print("Identifying differentially expressed features")
diff_expr_data <- get_diff_expr_data(expression_data, metadata, args)

print("Normalizing read count data")
norm_counts_data <- get_norm_counts_data(diff_expr_data$raw, metadata, args)     # includes all genes even those that have NA in diff_expr_data$res

export_plots(diff_expr_data, norm_counts_data, metadata, args)

print("Exporting differentially expressed features")
export_data(
    diff_expr_data$res,                                                          # this is not filtered differentially expressed features
    location=paste(args$output, "diff_expr_features.tsv", sep="_"),
    digits=5
)

print(
    paste(
        "Filtering normalized read counts matrix to include",
        "only differentially expressed features with padj <= ", args$padj
    )
)

row_metadata <- diff_expr_data$res %>%
                remove_rownames() %>%
                column_to_rownames("feature") %>%
                dplyr::select(log2FoldChange, pvalue, padj)  %>%                 # we are interested only in these three columns
                filter(.$padj<=args$padj) %>%
                arrange(desc(log2FoldChange))

col_metadata <- metadata %>%
                mutate_at(colnames(.), as.vector)                                # need to convert to vector, because in our metadata everything was a factor

norm_counts_mat <- assay(norm_counts_data)[as.vector(rownames(row_metadata)), ]
print("Size of the normalized read counts matrix after filtering")
print(dim(norm_counts_mat))

if (!is.null(args$cluster)){
    if (args$cluster == "column" || args$cluster == "both") {
        print("Clustering filtered read counts by columns")
        clustered_data = get_clustered_data(
            expression_data=norm_counts_mat,
            center=NULL,                                            # centering doesn't influence on the samples order
            dist="euclid",                                          # "euclidean distance is often a good choice for clustering arrays"
            transpose=TRUE
        )
        # norm_counts_mat <- clustered_data$expression              # no need to update norm_counts_mat as it should be the same because center=NULL
        col_metadata <- col_metadata[clustered_data$order, ]        # reordering samples order based on the HOPACH clustering resutls
        print("Reordered samples")
        print(col_metadata)
    }
    if (args$cluster == "row" || args$cluster == "both") {
        print("Clustering filtered normalized read counts by rows")
        clustered_data = get_clustered_data(
            expression_data=norm_counts_mat,
            center=if(args$center) "mean" else NULL,                # about centering normalized data https://www.biostars.org/p/387863/
            dist="cosangle",                                        # "cosangle distance metric is often a good choice for clustering genes"
            transpose=FALSE
        )
        norm_counts_mat <- clustered_data$expression                # can be different because of centering by rows mean
        row_metadata <- row_metadata[clustered_data$order, ]        # reordering features order based on the HOPACH clustering results
        print("Reordered features")
        print(head(col_metadata))
    }
}

print("Exporting normalized read counts to GCT format")
export_gct(
    counts_mat=norm_counts_mat,
    row_metadata=row_metadata,                                      # includes features as row names
    col_metadata=col_metadata,                                      # includes samples as row names
    location=paste(args$output, "_norm_read_counts.gct", sep="")
)