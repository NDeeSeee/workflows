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
suppressMessages(library(tidyverse))
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
        Idents(seurat_data) <- "new.ident"
        identity_data <- unique(as.vector(as.character(Idents(seurat_data))))
        if ( (nrow(extra_metadata) == length(identity_data)) && all(is.element(identity_data, extra_metadata$library_id)) ){
            print(paste("Extra metadata is successfully loaded from ", args$condition))
            print(extra_metadata)
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
        print("Extra metadata was not provided. Skipping.")
        return (seurat_data)
    }
}


get_filtered_data <- function(seurat_data, args){
    if (!is.null(args$groupby) && !is.null(args$select)){
        print(
            paste0(
                "Include only '",
                paste(args$select, collapse=" "), "' values from the '",
                args$groupby, "' metadata column, and only '",
                args$first, "' and '", args$second, "' values from the '",
                args$splitby, "' column."
            )
        )
        Idents(seurat_data) <- args$groupby
        seurat_data <- subset(seurat_data, idents=args$select)
        Idents(seurat_data) <- args$splitby
        seurat_data <- subset(seurat_data, idents=c(args$first, args$second))
        Idents(seurat_data) <- "new.ident"
    }
    return (seurat_data)
}


get_aggr_seurat_data <- function(seurat_data, group_by, features=NULL, assay="RNA", slot="counts"){
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay         # need to explicitely set RNA as scaling is run only on active assay
    Idents(seurat_data) <- group_by
    aggr_expr_data <- AggregateExpression(
        seurat_data,
        assays=assay,                          # need only RNA assay
        slot=slot,                             # for slot="counts" no exponentiation is performed prior to aggregating
        features=features,                     # if NULL use all features in assay
        return.seurat=TRUE,                    # summed values are saved in "counts", log-normalized - in "data", and scaled - in "scale.data"
        verbose=FALSE
    )
    Idents(seurat_data) <- "new.ident"
    DefaultAssay(seurat_data) <- backup_assay
    return (aggr_expr_data)
}


get_coldata <- function(seurat_data, args){
    col_data <- seurat_data@meta.data %>%
                dplyr::select(new.ident, all_of(args$splitby), any_of(args$batchby)) %>%     # use any_of(args$batchby) because it can be NULL
                distinct() %>%
                arrange(new.ident) %>%                                                       # no particular reasons to sort by rows 
                remove_rownames() %>%
                column_to_rownames("new.ident")
    col_data[[args$splitby]] <- as.factor(col_data[[args$splitby]])                          # DEseq likes factors
    if (!is.null(args$batchby)){
        col_data[[args$batchby]] <- as.factor(col_data[[args$batchby]])                      # DEseq likes factors
    }
    return (col_data)
}


get_diff_expr_data <- function(seurat_data, args, min_diff_pct=-Inf){
    all_features <- as.vector(as.character(rownames(seurat_data)))
    selected_features <- all_features[!all_features %in% args$exgenes]
    excluded_features <- all_features[all_features %in% args$exgenes]   # not used elsewhere
    if (length(excluded_features) > 0){
        print(
            paste(
                "Following genes will be excluded from the differential expression analysis:",
                paste(excluded_features, collapse=", ")
            )
        )
    }
    if (args$pseudo){
        print(
            paste(
                "Aggregating gene expression of the cells from the same datasets into",
                "the pseudobulk RNA-Seq samples. Using DESeq2 with sum of not normalized counts."
            )
        )
        aggr_seurat_data <- get_aggr_seurat_data(
            seurat_data,
            group_by="new.ident",                   # aggregating by dataset
            features=selected_features,
            assay="RNA",
            slot="counts"
        )
        counts <- GetAssayData(
            aggr_seurat_data,
            assay="RNA",                            # set assay in case GetAssayData won't take the default value
            slot="counts"                           # we need not normalized counts
        )
        col_data <- get_coldata(seurat_data, args)
        design_formula <- as.formula(
            paste0(
                "~",
                args$splitby,
                ifelse(
                    is.null(args$batchby),
                    "",
                    paste0("+", args$batchby)
                )
            )
        )
        deseq_data <- DESeqDataSetFromMatrix(
            countData=counts,
            colData=col_data,
            design=design_formula
        )
        if (args$lrt && !is.null(args$batchby)){
            print("Using LRT test to calculate p-values")           # see this post for details https://support.bioconductor.org/p/95493/#95572
            deseq_data <- DESeq(
                deseq_data,
                test="LRT",
                reduced=as.formula(paste0("~", args$splitby)),
                quiet=TRUE,
                parallel=TRUE
            )
        } else {
            print("Using Wald test to calculate p-values")
            deseq_data <- DESeq(
                deseq_data,
                quiet=TRUE,
                parallel=TRUE
            )
        }
        print(resultsNames(deseq_data))
        diff_expr_genes <- results(
            deseq_data,
            contrast=c(args$splitby, args$second, args$first),
            alpha=args$maxpvadj,                                 # recommended to set to our FDR threshold https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
            parallel=TRUE
        )
        print(mcols(diff_expr_genes))
        diff_expr_genes$log2FoldChange[is.na(diff_expr_genes$log2FoldChange)] <- 0.001         # we can't set it to 0, because it's log2
        diff_expr_genes$pvalue[is.na(diff_expr_genes$pvalue)] <- 1
        diff_expr_genes$padj[is.na(diff_expr_genes$padj)] <- 1
        diff_expr_genes <- as.data.frame(diff_expr_genes) %>%
                           dplyr::rename("avg_log2FC"="log2FoldChange") %>%
                           dplyr::rename("p_val"="pvalue") %>%
                           dplyr::rename("p_val_adj"="padj") %>%
                           rownames_to_column(var="gene") %>%
                           add_column(pct.1=0) %>%                                             # add pct.1 for consistency with Seurat's FindMarkers output
                           add_column(pct.2=0) %>%                                             # add pct.2 for consistency with Seurat's FindMarkers output
                           dplyr::select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
        return (list(diff_expr_genes=diff_expr_genes, data=deseq_data))
    } else {
        print(paste0(
            "Using Seurat FindMarkers functions with '", args$testuse, "' test for normalized counts")
        )
        Idents(seurat_data) <- args$splitby
        diff_expr_genes <- FindMarkers(
            seurat_data,
            assay="RNA",
            slot="data",                               # using normalized counts
            ident.1=args$second,
            ident.2=args$first,
            features=selected_features,
            logfc.threshold=args$minlogfc,
            min.pct=args$minpct,
            test.use=args$testuse,
            latent.vars=args$batchby,
            min.diff.pct=min_diff_pct,
            base=2,                                    # to make sure we use log2 scale
            verbose=TRUE
        ) %>% rownames_to_column(var="gene")
        Idents(seurat_data) <- "new.ident"
        return (list(diff_expr_genes=diff_expr_genes, data=seurat_data))                       # return here seurat_data just for consistency with output when we run DESeq2
    }
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


export_counts_plot <- function(data, col_data, from_deseq, features, rootname, splitby, batchby, x_label, y_label, legend_title, plot_title, combine_guides=NULL, palette="Paired", alpha=1, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plots = list()
            if(from_deseq){
                normalized_counts <- t(as.data.frame(assay(rlog(data, blind=FALSE))))    # transposed rlog normalized counts
                normalized_counts <- normalized_counts[, features, drop=FALSE]           # to include only features we need
                # if (splitby != batchby){                                               # remove batch effect
                #     normalized_counts <- vst(data, blind=FALSE)

                #     export_pca_plot(normalized_counts, "pca", c("condition"))
                #     mat <- limma::removeBatchEffect(assay(normalized_counts), batch=normalized_counts[[batchby]])
                #     assay(normalized_counts) <- mat
                #     export_pca_plot(normalized_counts, "pca_corrected", c("condition", "day"))
                #     normalized_counts <- t(as.data.frame(mat))
                # }
            } else {
                aggr_seurat_data <- get_aggr_seurat_data(
                    data,
                    group_by="new.ident",
                    features=features,
                    assay="RNA",
                    slot="counts"
                )
                normalized_counts <- t(                         # transposed log normalized counts
                    GetAssayData(
                        aggr_seurat_data,
                        assay="RNA",                            # set assay in case GetAssayData won't take the default value
                        slot="data"                             # we need normalized counts
                    )
                )
            }
            normalized_counts <- merge(normalized_counts, col_data, by=0) %>%      # merge by row names
                                 remove_rownames() %>%
                                 column_to_rownames("Row.names")
            for (i in 1:length(features)){
                current_feature <- features[i]
                counts <- normalized_counts[, c(current_feature, colnames(col_data)), drop=FALSE] %>%
                          dplyr::rename("count"=all_of(current_feature))
                plots[[i]] <- ggplot(
                                  counts,
                                  aes_string(x=batchby, y="count", color=splitby, group=splitby)
                              ) +
                              geom_point(alpha=alpha) +
                              stat_summary(fun=mean, geom="line") +
                              ggtitle(current_feature) +
                              theme_gray() +
                              xlab(x_label) +
                              ylab(y_label)
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


export_dim_plot <- function(data, rootname, reduction, plot_title, legend_title, perc_legend_title, split_by=NULL, group_by=NULL, perc_split_by=NULL, perc_group_by=NULL, label=FALSE, palette=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- DimPlot(
                        data,
                        reduction=reduction,
                        split.by=split_by,
                        group.by=group_by,
                        label=label
                    ) +
                    theme_gray() +
                    ggtitle(plot_title) +
                    guides(color=guide_legend(legend_title, override.aes=list(size=3)))

            if (!is.null(palette)){ plot <- plot + scale_color_brewer(palette=palette) }

            if(!is.null(perc_split_by) && !is.null(perc_group_by)){
                width <- 2 * width
                perc_data <- data@meta.data %>%
                             dplyr::group_by(across(all_of(perc_split_by)), across(all_of(perc_group_by))) %>%
                             summarise(counts=n(), .groups="drop") %>%               # drop the perc_group_by grouping level so we can get only groups defined by perc_split_by
                             mutate(freq=counts/sum(counts)*100) %>%                 # sum is taken for the group defined by perc_split_by
                             ungroup() %>%
                             complete(
                                 (!!as.symbol(perc_split_by)), (!!as.symbol(perc_group_by)),
                                 fill=list(counts=0, freq=0)
                             ) %>%
                             mutate(caption=paste0(counts, " (", round(freq, digits=2), "%)")) %>%
                             arrange(all_of(perc_split_by), all_of(perc_group_by))        # sort for consistency
                perc_data[[perc_group_by]] <- as.factor(perc_data[[perc_group_by]])
                label_data <- data@meta.data %>%
                              dplyr::group_by(across(all_of(perc_split_by))) %>%
                              summarise(counts=n(), .groups="drop_last") %>%              # drops all grouping as we have only one level
                              arrange(all_of(perc_split_by))                              # sort for consistency
                perc_plot <- ggplot(perc_data, aes_string(x=perc_split_by, y="freq", fill=perc_group_by)) +
                             geom_col(position="dodge", width=0.98, linetype="solid", color="black", show.legend=TRUE) +
                             geom_label(
                                aes(label=caption),
                                position=position_dodge(0.98),
                                color="white", size=4, show.legend=FALSE
                             ) +
                             guides(fill=guide_legend(perc_legend_title)) +
                             xlab("") +
                             ylab("Cells percentage") +
                             theme_gray() +
                             geom_label_repel(
                                label_data, mapping=aes(y=-Inf, label=counts),
                                color="black", fill="white", segment.colour=NA,
                                direction="y", size=4, show.legend=FALSE
                             ) +
                             RotatedAxis()
                plot <- plot + perc_plot
            }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export dim plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export dim plot to ", rootname, ".(png/pdf)", sep=""))
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
    parser <- ArgumentParser(
        description="Differential expression analysis for a subset of cells between two selected conditions"
    )
    parser$add_argument(
        "--rds",
        help=paste(
            "Path to the RDS file to load Seurat object from. RDS file can be produced by run_seurat.R script."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--condition",
        help=paste(
            "Path to the TSV/CSV file to optionally extend Seurat object metadata. First",
            "column 'library_id' should include all unique values from the 'new.ident'",
            "column of the loaded from --rds Seurat object metadata. All other columns will",
            "be added to the Seurat object metadata. If any of the provided in this file",
            "columns were already present in the Seurat object metadata, they will be",
            "overwritten.",
            "Default: no metadata columns will be added or overwritten"
        ),
        type="character"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from the Seurat object metadata to split cells into two groups",
            "to run --second vs --first differential expression analysis. May include",
            "columns from the metadata fields added with --condition."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--first",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby to define the",
            "first group of cells or pseudobulk RNA-Seq samples (when using --pseudo)."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--second",
        help=paste(
            "Value from the Seurat object metadata column set with --splitby to define the",
            "the second group of cells or pseudobulk RNA-Seq samples (when using --pseudo)"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--batchby",
        help=paste(
            "Column from the Seurat object metadata to define the variable that should",
            "be modelled as a batch effect when running differential expression analysis.",
            "Applied only when --testuse is one of 'LR', 'negbinom', 'poisson', or 'MAST',",
            "or when using --pseudo. May include columns from the metadata fields added",
            "with --condition. Values selected from the column set with --batchby should",
            "establish 1:1 relation with the 'new.ident' column of the Seurat object loaded",
            "from --rds.",
            "Default: do not model batch effect."
        ),
        type="character"
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Column from the Seurat object metadata to group cells for optional",
            "subsetting (for example, subset to the specific cluster or predicted",
            "cell type). May include columns from the metadata fields added with",
            "--condition."
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
            "change expression values. Ignored with --genes.",
            "Default: 10"
        ),
        type="integer", default=10
    )
    parser$add_argument(
        "--minlogfc",
        help=paste(
            "Include only those genes that on average have the absolute value of log2",
            "fold change expression difference not lower than this value. Increasing",
            "--minlogfc speeds up calculations, but can cause missing weaker signals.",
            "Ignored with --pseudo.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--minpct",
        help=paste(
            "Include only those genes that are detected in not lower than this fraction of cells",
            "in either of the two tested groups. Increasing --minpct speeds up calculations by not",
            "testing genes that are very infrequently expressed. Ignored with --pseudo.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--maxpvadj",
        help=paste(
            "Include only those genes for which adjusted P-val is not bigger that this value.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--testuse",
        help=paste(
            "Statistical test to use for differential gene expression analysis.",
            "Ignored with --pseudo.",
            "Default: wilcox"
        ),
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2"),
        type="character", default="wilcox"
    )
    parser$add_argument(
        "--pseudo",
        help=paste(
            "Aggregate gene expression of the cells from the same dataset into a pseudobulk",
            "RNA-Seq sample before running differential expression analysis with DESeq2.",
            "The following parameters will be ignored: --testuse, --minpct, --minlogfc.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--lrt",
        help=paste(
            "Use LRT instead of the pair-wise Wald test. Shows any differences across the variable",
            "set with --batchby whith the log2 fold changes calculated as the average expression",
            "changes due to criteria set with --splitby. Ignored when --pseudo or --batchby",
            "parameters are not provided.",
            "Default: use Wald test"
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

print(paste0("Loading Seurat data from '", args$rds, "'"))
seurat_data <- readRDS(args$rds)

print("Extending Seurat object with the extra metadata fields")
seurat_data <- extend_metadata(seurat_data, args)

print("Filtering Seurat object to include only required cells")
seurat_data <- get_filtered_data(seurat_data, args)

export_dim_plot(
    data=seurat_data,
    reduction="umap",
    plot_title=paste0("Split by '", args$splitby, "' clustered UMAP of filtered datasets"),
    legend_title=ifelse(
        is.null(args$groupby),
        args$splitby,
        args$groupby
    ),
    split_by=args$splitby,
    group_by=ifelse(
        is.null(args$groupby),
        args$splitby,
        args$groupby
    ),
    perc_split_by=args$splitby,
    perc_group_by=if(is.null(args$batchby)) args$groupby else args$batchby,                           # need to use this trick because ifelse doesn't return NULL
    perc_legend_title=if(is.null(args$batchby)) args$groupby else args$batchby,    # need to use this trick because ifelse doesn't return NULL
    label=ifelse(is.null(args$groupby), FALSE, TRUE),
    rootname=paste(args$output, "umap", sep="_"),
    pdf=args$pdf
)

print(
    paste0(
        "Running differential expression analysis for '", args$second, "' vs '", args$first,
        "' for cells or pseudobulk RNA-Seq samples split by '", args$splitby, "'"
    )
)

diff_expr_data <- get_diff_expr_data(seurat_data, args)

diff_expr_genes <- diff_expr_data$diff_expr_genes %>%
                   filter(.$p_val_adj<=args$maxpvadj) %>%
                   arrange(desc(avg_log2FC))

max_log2FC <- max(diff_expr_genes$avg_log2FC[is.finite(diff_expr_genes$avg_log2FC)])
min_log2FC <- min(diff_expr_genes$avg_log2FC[is.finite(diff_expr_genes$avg_log2FC)])
diff_expr_genes <- diff_expr_genes %>%                                                     # replace all possible Inf and -Inf with numbers
                   mutate("avg_log2FC"=ifelse(
                       .$avg_log2FC==Inf, max_log2FC+0.05*abs(max_log2FC), ifelse(
                           .$avg_log2FC==-Inf, min_log2FC-0.05*abs(min_log2FC), .$avg_log2FC)
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

export_volcano_plot(
    data=diff_expr_genes,
    rootname=paste(args$output, "diff_expr_genes", sep="_"),
    x_axis="avg_log2FC",
    y_axis="p_val_adj",
    x_cutoff=ifelse(args$pseudo, 0, args$minlogfc),
    y_cutoff=args$maxpvadj,
    x_label=bquote(~Log[2] ~ "FC"),
    y_label=bquote(~-Log[10] ~ italic(Padj)),
    plot_title="Differentially expressed genes",
    plot_subtitle=paste0(
        args$second, " vs ", args$first, " for ",
        ifelse(
            args$pseudo,
            "pseudobulk RNA-Seq samples",
            "cells"
        ),
        " split by ",
        args$splitby,
        " (",
        ifelse(
            args$pseudo,
            "",
            paste0("log2FC >= ", args$minlogfc, ", ")
        ),
        "Padj <= ", args$maxpvadj, ")",
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
            paste0("from '", paste(args$select, collapse=" "), "' subset(s) of '", args$groupby, "' metadata column"),
            ""
        )
    ),
    genes=highlight_genes,
    pdf=args$pdf
)

export_counts_plot(
    data=diff_expr_data$data,
    col_data=get_coldata(seurat_data, args),
    from_deseq=args$pseudo,                                               # when we use --pseudo the diff_expr_data$data comes from DESeq2
    features=highlight_genes,
    rootname=paste(args$output, "counts", sep="_"),
    splitby=args$splitby,
    batchby=ifelse(is.null(args$batchby), args$splitby, args$batchby),
    x_label=ifelse(is.null(args$batchby), args$splitby, args$batchby),
    y_label="Counts",
    legend_title=deparse(substitute(args$splitby)),
    plot_title=paste0(
        "Aggregated ",
        ifelse(args$pseudo, "rlog-", "log-"),
        "normalized counts"
    ),
    combine_guides="collect",
    pdf=args$pdf
)

print("Exporting differentially expressed genes")
export_data(diff_expr_genes, paste(args$output, "diff_expr_genes.tsv", sep="_"))