#!/usr/bin/env Rscript
options(warn = -1)
options("width" = 300)
options(
  error = function() {
    traceback(3)
    quit(
      save = "no",
      status = 1,
      runLast = FALSE
    )
  }
)

suppressMessages(library(argparse))
suppressMessages(library(BiocParallel))
suppressMessages(library(pheatmap))
suppressMessages(library(Glimma))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(hopach))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(DESeq2))
suppressMessages(library(limma)) # For removeBatchEffect
suppressMessages(library(cmapR))

mutate <- dplyr::mutate
filter <- dplyr::filter
group_by <- dplyr::group_by
slice <- dplyr::slice
rename <- dplyr::rename
select <- dplyr::select
arrange <- dplyr::arrange
distinct <- dplyr::distinct
`%>%` <- magrittr::`%>%`
`%in%` <- base::`%in%`
`%/%` <- base::`%/%`
`%%` <- base::`%%`

##########################################################################################
#
# v0.0.5 LRT Step 2
#
# Changes:
# - Applied 'limmaremovebatcheffect' batch correction in Step 2 if specified.
# - Used batch information from Step 1.
# - Ensured necessary data is passed between steps.
#
################################################################################################
# v0.0.3 - Modified to save all possible information into RDS files and handle batch correction
#
# Changes:
# - Added batch correction options (`CombatSeq` and `limmaRemoveBatchEffect`).
# - Modified design formula based on batch correction method.
# - Saved all necessary data (dds, contrasts, expression data, metadata) into RDS files.
#
##########################################################################################
#
# v0.0.2 LRT Step 2
#
# Changes:
# - Modified to read RDS files generated from the first script.
# - Inputs adjusted to accept RDS files and contrast indices.
# - Removed redundant code; now uses precomputed data.
#
##########################################################################################
#
# v0.0.1 LRT Step 2
#
# - Implemented **contrast generation** for complex experimental designs:
#   - Supports **main effects** (e.g., factor1 levels) and **interaction effects** (factor1 vs factor2).
#   - Automatically generates contrasts based on user-defined factors and levels.
#   - Designed for comparing multiple conditions and extracting meaningful differential gene expression patterns.
#
# - Added **log2 fold-change threshold** (default: 0.59, ~1.5 fold change) for filtering significant genes.
#
# - Integrated **batch correction** with two options:
#   - `combatseq`: Corrects for batch effects prior to differential analysis.
#   - `limmaremovebatcheffect`: Removes batch effects post-analysis using limma.
#
# - Visual outputs include:
#   - **Heatmaps** displaying top 30 most variable genes.
#   - **MDS plots** for visualizing sample distances.
#   - **MA plots** for visualizing differential expression.
#
# - **Clustering**: Supports HOPACH clustering with min-max scaling (row/column/both).
#
# - Exported **GCT files** compatible with GSEA, for both filtered and unfiltered data.
#
# - Multi-threading support through BiocParallel for faster processing.
#
# - Handles input in CSV/TSV formats and outputs standardized files.
#
##########################################################################################
#
# All input CSV/TSV files should have the following header (case-sensitive)
# <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
# <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV
#
# Format of the input files is identified based on file's extension
# *.csv - CSV
# *.tsv - TSV
# Otherwise used CSV by default
#
# The output file's rows order corresponds to the rows order of the first CSV/TSV file in
# the untreated group. Output is always saved in TSV format
#
# Output file includes only intersected rows from all input files. Intersected by
# RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand
#
# DESeq/DESeq2 always compares untreated_vs_treated groups (condition-1-vs-condition-2)
#
# Additionally we calculate -LOG10(pval) and -LOG10(padj)
#
# Use -un and -tn to set custom names for treated and untreated conditions
#
# Use -ua and -ta to set aliases for input expression files. Should be unique
# Exports GCT files to be used by GSEA. GCT files is always with uppercase GeneId
#
##########################################################################################


READ_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")

get_args <- function() {
  parser <- ArgumentParser(description = "Run DESeq2 analysis using contrasts from previous LRT step")
  parser$add_argument(
    "-e",
    "--expression_data_rds",
    help = "RDS file containing the expression data from step 1",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-c",
    "--contrasts_rds",
    help = "RDS file containing the contrasts list from step 1",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-w",
    "--dsq_wald_rds",
    help = "RDS file containing the DESeq2 object from the Wald test in step 1",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-m",
    "--metadata_rds",
    help = "RDS file containing the metadata from step 1",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-bcm",
    "--batch_correction_method_rds",
    help = "RDS file containing the batch correction method used in step 1",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-ci",
    "--contrast_indices",
    help = "Comma-separated list of integers representing contrast indices (e.g., 1,2,3)",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "--fdr",
    help = paste(
      "In the exploratory visualization part of the analysis output only features",
      "with adjusted p-value (FDR) not bigger than this value. Also the significance",
      "cutoff used for optimizing the independent filtering. Default: 0.1."
    ),
    type = "double",
    default = 0.1
  )
  parser$add_argument(
    "--lfcthreshold",
    help = paste(
      "Log2 fold change threshold for determining significant differential expression.",
      "Genes with absolute log2 fold change greater than this threshold will be considered.",
      "Default: 0.59 (about 1.5 fold change)"
    ),
    type = "double",
    default = 0.59
  )
  parser$add_argument(
    "--use_lfc_thresh",
    help = paste(
      "Flag to indicate whether to use lfcthreshold as the null hypothesis value in the results function call.",
      "If TRUE, lfcthreshold is used in the hypothesis test (i.e., genes are tested against this threshold).",
      "If FALSE, the null hypothesis is set to 0, and lfcthreshold is used only as a downstream filter.",
      "Default: FALSE"
    ),
    action = "store_true"
  )
  parser$add_argument(
    "--regulation",
    help = paste(
      "Direction of differential expression comparison. β is the log2 fold change.",
      "'both' for both up and downregulated genes (|β| > lfcThreshold); ",
      "'up' for upregulated genes (β > lfcThreshold); ",
      "'down' for downregulated genes (β < -lfcThreshold). ",
      "Default: both"
    ),
    type = "character",
    choices = c("both", "up", "down"),
    default = "both"
  )
  parser$add_argument(
    "-o",
    "--output",
    help = "Output prefix. Default: deseq",
    type = "character",
    default = "./deseq"
  )
  parser$add_argument(
    "-p",
    "--threads",
    help = "Threads",
    type = "integer",
    default = 1
  )
  parser$add_argument(
    "--test_mode",
    help = "Run for test, only first 100 rows",
    action = "store_true",
    default = FALSE
  )
  args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
  return(args)
}

# Load arguments
args <- get_args()

# Set threads
register(MulticoreParam(args$threads))

# Function to log messages with timestamps
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Load RDS files with detailed logging
log_message(paste("Loading expression data from", args$expression_data_rds))
expression_data_df <- readRDS(args$expression_data_rds)
log_message("Expression Data Loaded:")
print(head(expression_data_df))
log_message("Structure of Expression Data:")
glimpse(expression_data_df)

log_message(paste("Loading contrasts from", args$contrasts_rds))
all_contrasts <- readRDS(args$contrasts_rds)
log_message("Contrasts Loaded:")
print(all_contrasts)
log_message("Structure of Contrasts:")
glimpse(all_contrasts)

log_message(paste("Loading DESeq2 object from", args$dsq_wald_rds))
dds <- readRDS(args$dsq_wald_rds)
log_message("DESeq2 Object Loaded:")
print(dds)
log_message("Sample Names in DESeq2 Object:")
print(colnames(dds))

log_message(paste("Loading metadata from", args$metadata_rds))
metadata_df <- readRDS(args$metadata_rds)
log_message("Metadata Loaded:")
print(head(metadata_df))
log_message("Structure of Metadata:")
glimpse(metadata_df)

log_message(paste("Loading batch correction method from", args$batch_correction_method_rds))
batch_correction_method <- readRDS(args$batch_correction_method_rds)
log_message(paste("Batch correction method used:", batch_correction_method))

# Convert the comma-separated string into a vector of integers
contrast_vector <- as.integer(unlist(strsplit(args$contrast_indices, ",")))
log_message(paste("Contrast Indices Parsed:", paste(contrast_vector, collapse = ", ")))

# Ensure that the contrasts indices are valid
if (any(contrast_vector > length(all_contrasts$contrast_number))) {
  stop("One or more contrast indices are out of bounds.")
}

# Harmonize sample names
log_message("Starting harmonization of sample names...")

# Define READ_COL correctly based on your actual column naming
READ_COL <- "TotalReads" # Ensure this matches your read counts column pattern

# Extract read counts columns
read_counts_columns <- grep(
  paste(READ_COL, sep = ""),
  colnames(expression_data_df),
  value = TRUE,
  ignore.case = TRUE
)

log_message(paste("Columns matching READ_COL pattern ('", READ_COL, "'):", paste(read_counts_columns, collapse = ", ")))

read_counts_data_df <- expression_data_df[read_counts_columns]
log_message("Read Counts Data Extracted:")
print(head(read_counts_data_df))
log_message("Structure of Read Counts Data:")
glimpse(read_counts_data_df)

# Clean and standardize column names
log_message("Cleaning and standardizing column names...")
original_colnames <- colnames(read_counts_data_df)
log_message("Original Column Names:")
print(original_colnames)

cleaned_colnames <- sapply(original_colnames, function(s) {
  # Remove the last word if separated by space
  words <- unlist(strsplit(s, " ", fixed = TRUE))
  if (length(words) > 1) {
    paste(head(words, -1), collapse = "_")
  } else {
    s
  }
})

# Convert to lowercase and trim whitespace
cleaned_colnames <- tolower(cleaned_colnames)
cleaned_colnames <- trimws(cleaned_colnames)
cleaned_colnames <- gsub("[^[:alnum:]_]", "", cleaned_colnames)

log_message("Cleaned Column Names:")
print(cleaned_colnames)

# Assign cleaned column names back to the dataframe
colnames(read_counts_data_df) <- cleaned_colnames
log_message("Read Counts Data with Cleaned Column Names:")
print(head(read_counts_data_df))
log_message("Sample Names in Read Counts Data:")
print(colnames(read_counts_data_df))

log_message("Sample Names in DESeq2 Object After Cleaning:")
print(colnames(dds))
print("Checking for exact match between DESeq2 object and read counts data...")
print(dds)

# Check for exact match (regardless of order)
if (!setequal(colnames(dds), colnames(read_counts_data_df))) {
  # Identify mismatched samples
  missing_in_read_counts <- setdiff(colnames(dds), colnames(read_counts_data_df))
  missing_in_deseq <- setdiff(colnames(read_counts_data_df), colnames(dds))

  if (length(missing_in_read_counts) > 0) {
    log_message("Samples in DESeq2 object but not in read counts data:")
    print(missing_in_read_counts)
  }

  if (length(missing_in_deseq) > 0) {
    log_message("Samples in read counts data but not in DESeq2 object:")
    print(missing_in_deseq)
  }

  stop("Sample names in DESeq2 object and read counts data do not match.")
} else {
  # Reorder read counts columns to match DESeq2 object
  read_counts_data_df <- read_counts_data_df[, colnames(dds)]
  log_message("Reordered read counts data to match DESeq2 object sample names.")

  # Verify the reordering
  if (!all(colnames(dds) == colnames(read_counts_data_df))) {
    stop("Reordering failed: Sample names still do not match.")
  } else {
    log_message("Read counts data successfully reordered to match DESeq2 object.")
  }
}

# Apply limma batch correction if specified
if (batch_correction_method == "limmaremovebatcheffect" && "batch" %in% colnames(metadata_df)) {
  print("Applying limma::removeBatchEffect for batch correction")
  normCounts <- counts(dds, normalized = TRUE)
  batch <- metadata_df$batch
  design_matrix <- model.matrix(~conditions, data = metadata_df)
  normCounts_corrected <- removeBatchEffect(normCounts, batch = batch, design = design_matrix)
} else {
  normCounts_corrected <- counts(dds, normalized = TRUE)
}

# Function to get DESeq2 results for a specific contrast
get_contrast_res <- function(contrast_row) {
  # Determine altHypothesis based on regulation
  altHypothesis <- if (args$regulation == "up") {
    "greater"
  } else if (args$regulation == "down") {
    "less"
  } else {
    "greaterAbs"
  }

  lfcThreshold <- if (args$use_lfc_thresh) args$lfcthreshold else 0

  # DESeq2 object is already computed; we can directly extract results
  res <- results(dds,
                 name = contrast_row$contrast,
                 alpha = args$fdr,
                 lfcThreshold = lfcThreshold,
                 independentFiltering = TRUE,
                 altHypothesis = altHypothesis
  )
  return(res)
}

# Function to export MDS plot
export_mds_html_plot <- function(norm_counts_data, location) {
  tryCatch(
    expr = {
      htmlwidgets::saveWidget(
        glimmaMDS(
          x = norm_counts_data,
          groups = as.data.frame(metadata_df),
          labels = rownames(metadata_df)
        ),
        file = location
      )
    },
    error = function(e) {
      print(paste0("Failed to export MDS plot to ", location, " with error - ", e))
    }
  )
}

# Function to export GCT files
export_gct <- function(counts_mat,
                       row_metadata,
                       col_metadata,
                       location) {
  tryCatch(
    expr = {
      row_metadata <- row_metadata %>%
        rownames_to_column("id") %>%
        mutate_at("id", as.vector)
      col_metadata <- col_metadata %>%
        rownames_to_column("id") %>%
        mutate_at("id", as.vector)
      gct_data <- new("GCT",
                      mat = counts_mat[row_metadata$id, col_metadata$id],
                      # to guarantee the order and number of row/columns
                      rdesc = row_metadata,
                      cdesc = col_metadata
      )
      cmapR::write_gct(
        ds = gct_data,
        ofile = location,
        appenddim = FALSE
      )
      print(paste("Exporting GCT data to", location, sep = " "))
    },
    error = function(e) {
      print(paste("Failed to export GCT data to", location, sep = " "))
    }
  )
}

# Function to export CLS files
export_cls <- function(categories, location) {
  tryCatch(
    expr = {
      output_stream <- file(location, "w")
      on.exit(close(output_stream), add = TRUE)
      cat(
        paste(length(categories),
              length(levels(categories)),
              "1",
              sep = "\t"
        ),
        paste("#", paste(
          unique(as.character(categories)),
          collapse = "\t"
        ), sep = "\t"),
        paste(paste(as.character(categories),
                    collapse = "\t"
        ), sep = "\t"),
        file = output_stream,
        sep = "\n"
      )
      print(paste("Exporting CLS data to", location,
                  sep = " "
      ))
    },
    error = function(e) {
      print(paste("Failed to export CLS data to ", location, "with error - ", e,
                  sep = ""
      ))
    }
  )
}

# Function to generate clusters
get_clustered_data <- function(expression_data, transpose = FALSE) {
  if (transpose) {
    print("Transposing expression data")
    expression_data <- t(expression_data)
  }

  expression_data <- apply(
    expression_data,
    1,
    FUN = function(x) {
      scale_min_max(x)
    }
  )

  print("Running HOPACH")
  hopach_results <- hopach(expression_data)

  if (transpose) {
    print("Transposing expression data")
    expression_data <- t(expression_data)
  }

  print("Parsing cluster labels")
  clusters <- as.data.frame(hopach_results$clustering$labels)
  colnames(clusters) <- "label"
  clusters <- cbind(clusters, "HCL" = outer(clusters$label, 10^c((nchar(
    trunc(clusters$label)
  )[1] - 1):0), function(a, b) {
    paste0("c", a %/% b %% 10)
  }))
  clusters <- clusters[, c(-1), drop = FALSE]

  return(list(
    order = as.vector(hopach_results$clustering$order),
    expression = expression_data,
    clusters = clusters
  ))
}

# Function for min-max scaling
scale_min_max <- function(x,
                          min_range = -2,
                          max_range = 2) {
  min_val <- min(x)
  max_val <- max(x)
  scaled_x <-
    (x - min_val) / (max_val - min_val) * (max_range - min_range) + min_range
  return(scaled_x)
}

# Function to export MA plot
export_ma_plot <- function(data,
                           rootname,
                           width = 800,
                           height = 800,
                           resolution = 72) {
  tryCatch(
    expr = {
      png(
        filename = paste(rootname, ".png", sep = ""),
        width = width,
        height = height,
        res = resolution
      )
      plotMA(data)
      dev.off()

      pdf(
        file = paste(rootname, ".pdf", sep = ""),
        width = round(width / resolution),
        height = round(height / resolution)
      )
      plotMA(data)
      dev.off()

      cat(paste("\nExport MA-plot to ", rootname, ".(png/pdf)", "\n",
                sep =
                  ""
      ))
    },
    error = function(e) {
      dev.off()
      cat(paste(
        "\nFailed to export MA-plot to ",
        rootname,
        ".(png/pdf)",
        "\n",
        sep = ""
      ))
    }
  )
}

# Function to export heatmap
export_heatmap <- function(mat_data,
                           column_data,
                           rootname,
                           width = 800,
                           height = 800,
                           resolution = 72) {
  tryCatch(
    expr = {
      png(
        filename = paste(rootname, ".png", sep = ""),
        width = width,
        height = height,
        res = resolution
      )
      pheatmap(
        mat = mat_data,
        main = "Top 30 genes from VST normalized read counts",
        annotation_col = column_data,
        cluster_rows = FALSE,
        show_rownames = TRUE,
        cluster_cols = FALSE
      )
      dev.off()

      pdf(
        file = paste(rootname, ".pdf", sep = ""),
        width = round(width / resolution),
        height = round(height / resolution)
      )
      pheatmap(
        mat = mat_data,
        main = "Top 30 genes from VST normalized read counts",
        annotation_col = column_data,
        cluster_rows = FALSE,
        show_rownames = TRUE,
        cluster_cols = FALSE
      )
      dev.off()

      cat(paste(
        "\nExport expression heatmap to ",
        rootname,
        ".(png/pdf)",
        "\n",
        sep = ""
      ))
    },
    error = function(e) {
      dev.off()
      cat(
        paste(
          "\nFailed to export expression heatmap to ",
          rootname,
          ".(png/pdf)",
          "\n",
          sep = ""
        )
      )
    }
  )
}

# Function to clean DESeq2 results (handle NA values)
clean_deseq_results <- function(DESeqRes) {
  DESeqRes$log2FoldChange[is.na(DESeqRes$log2FoldChange)] <- 0
  DESeqRes$pvalue[is.na(DESeqRes$pvalue)] <- 1
  DESeqRes$padj[is.na(DESeqRes$padj)] <- 1
  return(DESeqRes)
}

# Function to add metadata columns and create the final results data frame
add_metadata_to_results <- function(collected_isoforms, DESeqRes, read_count_cols, digits) {
  collected_isoforms <- data.frame(
    cbind(collected_isoforms[, !colnames(collected_isoforms) %in% read_count_cols], DESeqRes),
    check.names = F,
    check.rows = F
  )
  collected_isoforms[, "'-LOG10(pval)'"] <- format(-log(as.numeric(collected_isoforms$pvalue), 10), digits = digits)
  collected_isoforms[, "'-LOG10(padj)'"] <- format(-log(as.numeric(collected_isoforms$padj), 10), digits = digits)
  return(collected_isoforms)
}

# Function to export DESeq2 report
export_deseq_report <- function(collected_isoforms, output_prefix) {
  collected_isoforms_filename <- paste0(output_prefix, "_gene_exp_table.tsv")
  write.table(
    collected_isoforms,
    file = collected_isoforms_filename,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  print(paste("Export DESeq report to", collected_isoforms_filename))
}

# Function to export normalized counts and filtered counts to GCT format
export_gct_data <- function(normCounts, collected_isoforms, column_data, output_prefix) {
  tryCatch(
    expr = {
      # Prepare row metadata
      row_metadata <- collected_isoforms %>%
        dplyr::mutate(GeneId = toupper(GeneId)) %>% # Ensure GeneId is uppercase
        dplyr::distinct(GeneId, .keep_all = TRUE) %>% # Remove duplicates based on GeneId
        remove_rownames() %>%
        column_to_rownames("GeneId") %>% # Set GeneId as rownames
        dplyr::select(log2FoldChange, pvalue, padj) %>% # Select relevant columns
        arrange(desc(log2FoldChange)) # Arrange by log2FoldChange

      # Check for any remaining NAs
      if (any(is.na(rownames(row_metadata)))) {
        stop("There are still NA GeneIds after processing.")
      }

      # Prepare column metadata
      col_metadata <- column_data %>%
        mutate_at(colnames(.), as.vector)

      # Create GCT object
      gct_data <- new("GCT",
                      mat = normCounts,
                      rdesc = row_metadata,
                      cdesc = col_metadata
      )

      # Write GCT files
      cmapR::write_gct(
        ds = gct_data,
        ofile = paste0(output_prefix, "_counts_all.gct"),
        appenddim = FALSE
      )
      print(paste("Exporting GCT data to", paste0(output_prefix, "_counts_all.gct")))

      # Filter rows by padj <= FDR
      row_metadata_filtered <- row_metadata %>%
        dplyr::filter(padj <= args$fdr)

      # Check if any genes remain after filtering
      if (nrow(row_metadata_filtered) == 0) {
        warning(paste("No genes passed the FDR threshold of", args$fdr, "for", output_prefix))
      }

      # Subset normCounts to include only the filtered GeneIds
      filtered_normCounts <- normCounts[rownames(row_metadata_filtered),]

      # Create filtered GCT object
      gct_data_filtered <- new("GCT",
                               mat = filtered_normCounts,
                               rdesc = row_metadata_filtered,
                               cdesc = col_metadata
      )

      # Write filtered GCT file
      cmapR::write_gct(
        ds = gct_data_filtered,
        ofile = paste0(output_prefix, "_counts_filtered.gct"),
        appenddim = FALSE
      )
      print(paste("Exporting GCT data to", paste0(output_prefix, "_counts_filtered.gct")))
    },
    error = function(e) {
      print(paste("Failed to export GCT data to", output_prefix, "with error -", e$message))
    }
  )
}

# Function to cluster data and re-order based on clustering results
cluster_and_reorder <- function(normCounts, col_metadata, row_metadata, args) {
  if (!is.null(args$cluster)) {
    if (args$cluster == "column" || args$cluster == "both") {
      print("Clustering filtered read counts by columns")
      clustered_data <- get_clustered_data(expression_data = normCounts, transpose = TRUE)
      col_metadata <- cbind(col_metadata, clustered_data$clusters) # Add cluster labels
      col_metadata <- col_metadata[clustered_data$order,] # Reorder based on clustering results
      print("Reordered samples")
      print(col_metadata)
    }
    if (args$cluster == "row" || args$cluster == "both") {
      print("Clustering filtered normalized read counts by rows")
      clustered_data <- get_clustered_data(expression_data = normCounts, transpose = FALSE)
      normCounts <- clustered_data$expression # Adjust for row centering
      row_metadata <- cbind(row_metadata, clustered_data$clusters) # Add cluster labels
      row_metadata <- row_metadata[clustered_data$order,] # Reorder based on clustering results
      print("Reordered features")
      print(head(row_metadata))
    }
  }
  return(list(normCounts = normCounts, col_metadata = col_metadata, row_metadata = row_metadata))
}

# Function to export results (plots, reports, etc.)
export_contrast_results <- function(res, final_isoforms, contrast_name, dds, normCounts, output, args) {
  # Generate and export MA plot
  export_ma_plot(res, paste0(output, "_", contrast_name, "_ma_plot"))

  # Heatmap for the top 30 most expressed genes
  vst <- varianceStabilizingTransformation(dds, blind = FALSE)
  vsd <- assay(vst)
  rownames(vsd) <- rownames(normCounts)
  mat <- vsd[order(rowMeans(normCounts), decreasing = TRUE)[1:30],]
  export_heatmap(mat, as.data.frame(colData(dds)), paste0(output, "_", contrast_name, "_heatmap"))

  # Export the DESeq2 results to a TSV file
  results_filename <- paste0(output, "_", contrast_name, "_results.tsv")
  write.table(as.data.frame(res), file = results_filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  print(paste("Results for", contrast_name, "exported to", results_filename))

  # Export MDS plot
  export_mds_html_plot(normCounts, paste0(output, "_", contrast_name, "_mds_plot.html"))

  # Export GCT data
  export_gct_data(normCounts, final_isoforms, as.data.frame(colData(dds)), paste0(output, "_", contrast_name))
}

# Iterate through each index in contrast_vector and process the corresponding contrast
for (contrast_index in contrast_vector) {
  contrast_row <- all_contrasts[all_contrasts$contrast_number == contrast_index,]
  if (nrow(contrast_row) == 0) {
    print(paste("Contrast index", contrast_index, "not found in contrasts list. Skipping."))
    next
  }

  contrast_name <- paste0("contrast_", contrast_index, "_", contrast_row$contrast)
  print(paste("Processing contrast:", contrast_name))

  # Get DESeq2 results for the specific contrast
  res <- get_contrast_res(contrast_row)

  print("DESeq2 results obtained.")
  print("Summary:")
  print(summary(res))
  print("Head:")
  print(head(res))

  # Clean results
  DESeqRes <- clean_deseq_results(as.data.frame(res[, c(1, 2, 5, 6)]))

  print("DESeq2 results cleaned.")
  print(head(DESeqRes))

  expression_data_df <- expression_data_df %>%
    dplyr::mutate_at("GeneId", toupper) %>%
    dplyr::distinct(GeneId, .keep_all = TRUE)

  print("Expression data cleaned.")
  print(head(expression_data_df))

  # Add metadata to results
  final_isoforms <- add_metadata_to_results(expression_data_df, DESeqRes, read_counts_columns, digits = 4)

  print("Metadata added to results.")

  # Export DESeq2 report
  export_deseq_report(final_isoforms, paste0(args$output, "_", contrast_name))
  print("Exported DESeq2 report.")

  # **Removed** the direct call to export_gct_data here
  # print("Exported normalized counts to GCT format.")

  # Generate and export plots (including GCT export)
  export_contrast_results(res, final_isoforms, contrast_name, dds, normCounts_corrected, args$output, args)
}

print("DESeq2 analysis complete.")
