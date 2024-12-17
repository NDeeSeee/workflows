#!/usr/bin/env Rscript
options(warn = -1)
options("width" = 300)

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
suppressMessages(library(rlang)) # For traceback handling

options(rlang_backtrace_on_error = "full")

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
    "-dsq",
    "--dsq_obj_data",
    help = "RDS file containing contrassts and expression data from step 1",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-bcm",
    "--batchcorrection",
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
    "--cluster",
    help    = paste(
      "Hopach clustering method to be run on normalized read counts for the",
      "exploratory visualization part of the analysis. Default: do not run",
      "clustering"
    ),
    type    = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
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
    help   = "Run for test, only first 100 rows",
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
log_message(paste("Loading contrasts from", args$dsq_obj_data))
all_contrasts <- readRDS(args$dsq_obj_data)
log_message("Contrasts Loaded:")
print(all_contrasts)
log_message("Structure of Contrasts:")
glimpse(all_contrasts)

log_message(paste("Loading expression data from", args$dsq_obj_data))
expression_data_df <- all_contrasts$expression_data_df

all_contrasts <- purrr::list_modify(all_contrasts, expression_data_df = NULL)

log_message("Expression Data Loaded:")
print(head(expression_data_df))
log_message("Structure of Expression Data:")
glimpse(expression_data_df)

log_message(paste("Loading DESeq2 object from Contrasts"))
dds <- all_contrasts[[1]]$subset
log_message("DESeq2 Object Loaded:")
print(dds)
log_message("Sample Names in DESeq2 Object:")
print(colnames(dds))

log_message(paste("Loading metadata from DESeq2 object"))
metadata_df <- colData(dds)
log_message("Metadata Loaded:")
print(head(metadata_df))
log_message("Structure of Metadata:")
glimpse(metadata_df)

batch_correction_method <- args$batchcorrection
log_message(paste("Batch correction method used:", batch_correction_method))

# Convert the comma-separated string into a vector of integers
contrast_vector <- as.integer(unlist(strsplit(args$contrast_indices, ",")))
log_message(paste("Contrast Indices Parsed:", paste(contrast_vector, collapse = ", ")))

# Ensure that the contrasts indices are valid
if (any(contrast_vector > length(all_contrasts))) {
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
  # Case 3: After DESeq2, apply rlog transformation and remove batch effects using limma
  print("Applying rlog transformation and limma batch effect removal")
  rlog_transformed <- rlog(dds, blind = FALSE)
  rlog_counts <- assay(rlog_transformed)

  # Prepare design matrix without 'batch' for removeBatchEffect
  design_formula <- as.formula(paste0("~", str_remove(as.character(dds@design)[2], " \\+ batch")))

  design_matrix <- model.matrix(design_formula, data = metadata_df)

  # Apply removeBatchEffect
  normCounts <- limma::removeBatchEffect(rlog_counts, batch = metadata_df$batch, design = design_matrix)
} else {
  normCounts <- rlog(dds, blind = FALSE)
  normCounts <- assay(normCounts)
}

# Function to get DESeq2 results for a specific contrast
get_contrast_res <- function(contrast_row) {
  dds_subset <- contrast_row$subset

  # Print contrast details
  print(paste("DESeq object Subset for Contrast:", contrast_row$contrast))
  print(dds_subset)
  print(resultsNames(dds_subset))

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
  res <- results(dds_subset,
                 name                 = contrast_row$contrast,
                 alpha                = args$fdr,
                 lfcThreshold         = lfcThreshold,
                 independentFiltering = TRUE,
                 altHypothesis        = altHypothesis
  )

  print("DESeq2 subset results obtained.")
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

  # Apply scaling per row without changing dimensions
  expression_data <- t(apply(
    expression_data,
    1,
    FUN = function(x) {
      scale_min_max(x)
    }
  ))

  if (transpose) {
    print("Transposing expression data back")
    expression_data <- t(expression_data)
  }

  print("Running HOPACH")
  hopach_results <- hopach::hopach(expression_data, verbose = TRUE)

  print("Parsing cluster labels")
  clusters <- data.frame(label = hopach_results$clustering$labels)
  clusters <- clusters %>% mutate(HCL = paste0("c", label))

  return(list(
    order = hopach_results$clustering$order,
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

# Function to cluster data and re-order based on clustering results
cluster_and_reorder <- function(normCounts, col_metadata, row_metadata, args) {
  if (args$cluster != "none") { # Only cluster if not "none"
    if (args$cluster == "column" || args$cluster == "both") {
      print("Clustering filtered read counts by columns")
      clustered_data <- get_clustered_data(expression_data = normCounts, transpose = TRUE)
      col_metadata <- cbind(col_metadata, clustered_data$clusters)
      col_metadata <- col_metadata[clustered_data$order,]
      normCounts   <- normCounts[, clustered_data$order]
      print("Reordered samples")
      print(col_metadata)
    }
    if (args$cluster == "row" || args$cluster == "both") {
      print("Clustering filtered normalized read counts by rows")
      clustered_data <- get_clustered_data(expression_data = normCounts, transpose = FALSE)
      normCounts   <- clustered_data$expression
      row_metadata <- cbind(row_metadata, clustered_data$clusters)
      row_metadata <- row_metadata[clustered_data$order,]
      normCounts   <- normCounts[clustered_data$order,]
      print("Reordered features")
      print(head(row_metadata))
    }
  } else {
    print("Clustering skipped as per 'none' option.")
  }
  return(list(normCounts = normCounts, col_metadata = col_metadata, row_metadata = row_metadata))
}

# TODO: the following funcion should be updated to handle the any contrast data with modifiee LFC and padj column
#  names into appropriate contrast like contrast_1_LFC and so on and the filtering should be then based on either of
#  them

# Function to add metadata columns and create the final results data frame
add_metadata_to_results <- function(expression_data_df, deseq_result, read_count_cols, digits) {
  deseq_result <- deseq_result %>%
    as.data.frame() %>%
    select(baseMean, log2FoldChange, pvalue, padj) %>%
    # To handle NA values in pvalue and padj
    # mutate(
    #   `-LOG10(pval)` = -log10(as.numeric(pvalue)),
    #   `-LOG10(padj)` = -log10(as.numeric(padj))
    # )

    expression_data_df <- data.frame(
    cbind(expression_data_df[, !colnames(expression_data_df) %in% read_count_cols], deseq_result),
    check.names = F,
    check.rows = F
  )
  return(expression_data_df)
}

# Function to export DESeq2 report
export_deseq_report <- function(expression_data_df, output_prefix) {
  expression_data_df_filename <- paste0(output_prefix, "_gene_exp_table.tsv")
  write.table(
    expression_data_df,
    file = expression_data_df_filename,
    sep   = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  print(paste("Export DESeq report to", expression_data_df_filename))
}

# TODO: here filtering should be based on either of LFC or padj column for ANY of the contrast (like contrast_1_padj
#  <= fdr | contrast_2_padj <= fdr, etc.)
# Function to export normalized counts and filtered counts to GCT format
export_gct_data <- function(normCounts, row_metadata, col_metadata, output_prefix) {
  tryCatch(
    expr = {
      # Use the provided row_metadata directly
      # Prepare column metadata
      col_metadata <- col_metadata %>% mutate_all(as.vector)

      # Create GCT object
      gct_data <- new("GCT",
                      mat   = normCounts,
                      rdesc = row_metadata,
                      cdesc = col_metadata)

      # Write GCT files
      cmapR::write_gct(
        ds = gct_data,
        ofile = paste0(output_prefix, "_counts_all.gct"),
        appenddim = FALSE
      )
      print(paste("Exporting GCT data to", paste0(output_prefix, "_counts_all.gct")))

      # Filter rows by padj <= FDR if padj is available
      if ('padj' %in% colnames(row_metadata)) {
        row_metadata_filtered <- row_metadata %>% dplyr::filter(padj <= args$fdr)
      } else {
        row_metadata_filtered <- row_metadata
      }

      # Check if any genes remain after filtering
      if (nrow(row_metadata_filtered) == 0) {
        warning(paste("No genes passed the FDR threshold of", args$fdr, "for", output_prefix))
      }

      # Subset normCounts to include only the filtered rows
      filtered_normCounts <- normCounts[rownames(row_metadata_filtered),]

      # Create filtered GCT object
      gct_data_filtered <- new("GCT",
                               mat   = filtered_normCounts,
                               rdesc = row_metadata_filtered,
                               cdesc = col_metadata)

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

# Iterate through each index in contrast_vector and process the corresponding contrast
annotated_expression_combined_df <- data.frame()

for (contrast_index in contrast_vector) {
  contrast_selected <- all_contrasts[[contrast_index]]

  print("Selected contrast structure:")
  print(contrast_selected)

  contrast_name <- paste0("contrast_", contrast_index, "_", contrast_selected$contrast)
  print(paste("Processing contrast:", contrast_name))

  # Get DESeq2 results for the specific contrast
  deseq_result <- get_contrast_res(contrast_selected)

  print("DESeq2 results obtained.")
  print("Summary:")
  print(summary(deseq_result))
  print("Head:")
  print(head(deseq_result))

  expression_data_df <- expression_data_df %>%
    dplyr::mutate_at("GeneId", toupper) %>%
    dplyr::distinct(GeneId, .keep_all = TRUE)

  print("Expression data cleaned.")
  print(head(expression_data_df))

  # Add metadata to results
  annotated_expression_df <- add_metadata_to_results(expression_data_df, deseq_result, read_counts_columns, digits = 4)

  print("Metadata added to results.")

  # Export DESeq2 report
  export_deseq_report(annotated_expression_df, paste0(args$output, "_", contrast_name))
  print("Exported DESeq2 report.")

  # Generate and export plots (including GCT export)
  annotated_expression_combined_df <- bind_rows(annotated_expression_combined_df, annotated_expression_df)
}

export_mds_html_plot(normCounts, paste0(args$output, "_mds_plot.html"))

# After the loop, ensure that 'row_metadata' is correctly set
row_metadata <- annotated_expression_combined_df %>%
  dplyr::select(GeneId, log2FoldChange, pvalue, padj) %>%
  dplyr::distinct(GeneId, .keep_all = TRUE) %>%
  remove_rownames() %>%
  column_to_rownames("GeneId")

# TODO: here normCounts should be the filtered data
clustered_data <- cluster_and_reorder(normCounts, as.data.frame(colData(dds)), row_metadata, args)

# TODO: here normCounts should be the entire data
# Export GCT data
export_gct_data(clustered_data$normCounts, clustered_data$row_metadata, clustered_data$col_metadata, args$output)

print("DESeq2 analysis complete.")

