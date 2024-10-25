#!/usr/bin/env Rscript

# Load required libraries
suppressMessages({
  library(cmapR)
  library(argparse)
  library(tidyverse)  # Includes dplyr, purrr, tibble, magrittr, etc.
  library(stringr)
})

# Explicitly assign the pipe operator from magrittr to ensure consistency
`%>%` <- magrittr::`%>%`
`%in%` <- base::`%in%`

# Define command-line arguments using 'argparse'
parser <- ArgumentParser(description = "Merge multiple GCT files into one with per-contrast metadata.")

parser$add_argument("-i", "--input", required = TRUE, nargs = "+",
                    help = "List of GCT files to merge, each representing a different contrast.")

parser$add_argument("-o", "--output", default = "merged_counts.gct",
                    help = "Output filename [default: %(default)s].")

args <- parser$parse_args()

# Output the parsed arguments
cat("Parsed arguments:\n")
print(args)
flush.console()

# Function to extract contrast name from filename
extract_contrast_name <- function(filename) {
  # Example filename: deseq_output_contrast_11_counts_filtered.gct
  # Extract 'contrast_11' or any pattern after 'contrast_'
  match <- str_match(basename(filename), "contrast_(\\S+?)_counts_filtered\\.gct")
  if (is.na(match[1, 2])) {
    stop(paste("Cannot extract contrast name from filename:", filename))
  }
  return(paste0("contrast_", match[1, 2]))
}

# Initialize lists to store data
contrast_rdesc_list <- list()
gene_id_list <- list()
expression_df_list <- list()

# Read each GCT file and store its components
cat("Reading GCT files and extracting data...\n")
flush.console()

for (gct_file in args$input) {
  contrast_name <- extract_contrast_name(gct_file)
  cat("\nProcessing GCT file for", contrast_name, ":", gct_file, "\n")
  flush.console()

  # Read GCT file
  gct_data <- cmapR::parse_gctx(gct_file, matrix_only = FALSE)

  cat("Successfully read GCT file:", gct_file, "\n")
  flush.console()

  # Print sample names from the current GCT file
  current_sample_names <- colnames(gct_data@mat)
  cat("Sample names in current GCT file:\n")
  print(current_sample_names)
  flush.console()

  # Extract log2FoldChange and padj columns
  if (!all(c("log2FoldChange", "padj") %in% colnames(gct_data@rdesc))) {
    stop(paste("GCT file", gct_file, "does not contain required row metadata columns 'log2FoldChange' and 'padj'."))
  }

  # Extract log2FoldChange and padj, rename them with contrast name
  rdesc_contrast <- gct_data@rdesc %>%
    select(id, log2FoldChange, padj) %>%
    mutate(
      log2FoldChange = round(log2FoldChange, digits = 2),
      padj = if_else(
        abs(padj) < 0.001,
        sprintf("%.1e", padj),
        sprintf("%.3f", padj)
      )
    ) %>%
    rename(
      !!paste0("log2FoldChange_", contrast_name) := log2FoldChange,
      !!paste0("padj_", contrast_name) := padj
    )

  # Store the contrast-specific rdesc
  contrast_rdesc_list[[contrast_name]] <- rdesc_contrast

  # Store gene IDs
  gene_id_list[[contrast_name]] <- toupper(rdesc_contrast$id)

  # Extract expression data
  expr_df <- as.data.frame(gct_data@mat) %>%
    rownames_to_column("id") %>%
    mutate(id = toupper(id))

  # Store expression data
  expression_df_list[[contrast_name]] <- expr_df

  # Store column descriptors and sample names from the first GCT file
  if (!exists("merged_cdesc")) {
    merged_cdesc <- gct_data@cdesc
    # Remove 'id' column if present
    if ("id" %in% colnames(merged_cdesc)) {
      merged_cdesc <- merged_cdesc %>% select(-id)
    }
    # Ensure that row names of cdesc are sample names
    if (!all(rownames(merged_cdesc) == colnames(gct_data@mat))) {
      stop("Row names of column descriptors do not match sample names in the expression matrix.")
    }
    # Store sample names for validation
    sample_names <- colnames(gct_data@mat)
    cat("\nSample names from the first GCT file (used as reference):\n")
    print(sample_names)
    flush.console()
  } else {
    # Validate that sample names are consistent (using set equality)
    if (!setequal(colnames(gct_data@mat), sample_names)) {
      cat("\nMismatch in sample names detected.\n")
      cat("Sample names in reference:\n")
      print(sample_names)
      cat("Sample names in current GCT file:\n")
      print(colnames(gct_data@mat))
      stop(paste("Sample names do not match across GCT files. Mismatch found in file:", gct_file))
    }
  }
}

# Merge gene IDs from all contrasts
all_gene_ids <- unique(unlist(gene_id_list))

# Create a data frame with all gene IDs
merged_rdesc <- data.frame(id = all_gene_ids, stringsAsFactors = FALSE)

# Merge contrast-specific row metadata into merged_rdesc
cat("\nMerging contrast-specific row metadata...\n")
flush.console()

for (contrast_name in names(contrast_rdesc_list)) {
  rdesc_contrast <- contrast_rdesc_list[[contrast_name]]
  # Convert gene IDs to uppercase for consistency
  rdesc_contrast$id <- toupper(rdesc_contrast$id)
  # Merge into merged_rdesc
  merged_rdesc <- merged_rdesc %>%
    full_join(rdesc_contrast, by = "id")
}

# Remove 'pvalue' columns if they exist (since only 'padj' is required)
merged_rdesc <- merged_rdesc %>%
  select(-contains("pvalue"))

# Set 'id' as row names
merged_rdesc <- merged_rdesc %>%
  arrange(id) %>%
  remove_rownames() %>%
  column_to_rownames("id")

# Merge expression data from all GCT files
cat("\nMerging expression data from all GCT files...\n")
flush.console()

# Convert each expression data frame to long format
expression_long_list <- list()
for (contrast_name in names(expression_df_list)) {
  expr_df <- expression_df_list[[contrast_name]]
  expr_long <- expr_df %>%
    pivot_longer(
      cols = -id,
      names_to = "sample",
      values_to = "expression"
    )
  expression_long_list[[contrast_name]] <- expr_long
}

# Combine all expression data into one long data frame
expression_long_combined <- bind_rows(expression_long_list)

# Remove duplicates by taking the first non-NA expression value
expression_long_combined <- expression_long_combined %>%
  group_by(id, sample) %>%
  summarize(expression = first(na.omit(expression)), .groups = "drop")

# Pivot back to wide format
expression_df <- expression_long_combined %>%
  pivot_wider(
    names_from = sample,
    values_from = expression
  )

# Replace NAs with zeros
expression_df[is.na(expression_df)] <- 0

# Ensure that sample columns are consistent
expression_samples <- colnames(expression_df)[!colnames(expression_df) %in% c("id")]

# Print expression samples and reference sample names
cat("\nSample names in merged expression data:\n")
print(expression_samples)
cat("\nReference sample names from the first GCT file:\n")
print(sample_names)
flush.console()

# Check if sample names match (using set equality)
if (!setequal(expression_samples, sample_names)) {
  cat("\nMismatch in sample names after merging expression data.\n")
  cat("Sample names in expression data:\n")
  print(expression_samples)
  cat("Expected sample names:\n")
  print(sample_names)
  stop("Sample names in expression data do not match the expected sample names.")
}

# Reorder columns in expression_df to match sample_names
expression_df <- expression_df %>%
  select(id, all_of(sample_names))

# Convert back to matrix
expression_mat <- expression_df %>%
  remove_rownames() %>%
  column_to_rownames("id") %>%
  as.matrix()

# Ensure that row names of merged_rdesc and expression_mat match
merged_rdesc <- merged_rdesc[rownames(expression_mat),]

# Validate that the number of genes and their order match
if (!all(rownames(merged_rdesc) == rownames(expression_mat))) {
  stop("Mismatch between gene IDs in row metadata and expression matrix.")
}

# Prepare column metadata (merged_cdesc)
# Ensure that the number of samples matches between cdesc and expression matrix
num_samples <- ncol(expression_mat)
if (nrow(merged_cdesc) != num_samples) {
  stop(paste("Number of samples in cdesc (", nrow(merged_cdesc), ") does not match number of samples in the expression matrix (", num_samples, ")."))
}

# Ensure that sample names match between cdesc and expression matrix
if (!all(rownames(merged_cdesc) == colnames(expression_mat))) {
  cat("\nMismatch in sample names between column descriptors and expression matrix.\n")
  cat("Sample names in column descriptors:\n")
  print(rownames(merged_cdesc))
  cat("Sample names in expression matrix:\n")
  print(colnames(expression_mat))
  stop("Sample names in column descriptors do not match those in the expression matrix.")
}

# Prepare column metadata (ensure columns are vectors)
merged_cdesc <- merged_cdesc %>%
  mutate_all(as.vector)

# Create GCT object
cat("\nCreating GCT object...\n")
flush.console()

gct_obj <- methods::new("GCT",
                        mat = expression_mat,
                        rdesc = merged_rdesc,
                        cdesc = merged_cdesc,
                        rid = rownames(expression_mat),
                        cid = colnames(expression_mat),
                        src = args$output)

cat("Writing merged GCT file to", args$output, "...\n")
flush.console()

# Write GCT file
cmapR::write_gct(gct_obj, ofile = args$output, appenddim = FALSE)

cat("Merged GCT file successfully written to", args$output, "\n")
flush.console()