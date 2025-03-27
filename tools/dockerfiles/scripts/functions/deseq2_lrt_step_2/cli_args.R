#!/usr/bin/env Rscript
#
# Command-line argument handling for DESeq2 LRT Step 2
#

#' Define and parse command line arguments
#'
#' @return Parsed arguments list
#' @export
get_args <- function() {
  parser <- ArgumentParser(description = "Run DESeq2 analysis using contrasts from previous LRT step")
  parser$add_argument(
    "--dsq_obj_data",
    help = "RDS file containing contrasts and expression data from step 1",
    type = "character"
  )
  parser$add_argument(
    "--contrast_df",
    help = "TSV file containing contrasts data",
    type = "character"
  )
  parser$add_argument(
    "--batchcorrection",
    help = "Batch correction method to use: 'none', 'combatseq', or 'model'",
    type = "character",
    choices = c("none", "combatseq", "model"),
    default = "none"
  )
  parser$add_argument(
    "--contrast_indices",
    help = "Comma-separated list of integers representing contrast indices (e.g., 1,2,3)",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "--fdr",
    help = paste(
      "FDR cutoff for significance filtering. Default: 0.1."
    ),
    type = "double",
    default = 0.1
  )
  parser$add_argument(
    "--lfcthreshold",
    help = paste(
      "Log2 fold change threshold for determining significant differential expression.",
      "Default: 0.59 (about 1.5 fold change)"
    ),
    type = "double",
    default = 0.59
  )
  parser$add_argument(
    "--use_lfc_thresh",
    help = paste(
      "Use lfcthreshold as the null hypothesis value in the results function call.",
      "Default: FALSE"
    ),
    action = "store_true",
    default = FALSE
  )
  parser$add_argument(
    "--regulation",
    help = paste(
      "Direction of differential expression comparison: 'both', 'up', or 'down'.",
      "Default: both"
    ),
    type = "character",
    choices = c("both", "up", "down"),
    default = "both"
  )
  parser$add_argument(
    "--cluster",
    help = paste(
      "Hopach clustering method to be run on normalized read counts.",
      "Default: none"
    ),
    type = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
  )
  parser$add_argument(
    "--scaling_type",
    help = paste(
      "Type of scaling for expression data: 'minmax' or 'zscore'.",
      "Default: zscore"
    ),
    type = "character",
    choices = c("minmax", "zscore"),
    default = "zscore"
  )
  parser$add_argument(
    "--rowdist",
    help = paste(
      "Distance metric for HOPACH row clustering.",
      "Default: cosangle"
    ),
    type = "character",
    choices = c(
      "cosangle",
      "abscosangle",
      "euclid",
      "cor",
      "abscor"
    ),
    default = "cosangle"
  )
  parser$add_argument(
    "--columndist",
    help = paste(
      "Distance metric for HOPACH column clustering.",
      "Default: euclid"
    ),
    type = "character",
    choices = c(
      "cosangle",
      "abscosangle",
      "euclid",
      "cor",
      "abscor"
    ),
    default = "euclid"
  )
  parser$add_argument(
    "--k",
    help = "Number of levels for Hopach clustering (1-15). Default: 3.",
    type = "integer",
    default = 3
  )
  parser$add_argument(
    "--kmax",
    help = "Maximum number of clusters at each level (2-9). Default: 5.",
    type = "integer",
    default = 5
  )
  parser$add_argument(
    "--output",
    help = "Output prefix. Default: deseq-lrt-step-2",
    type = "character",
    default = "deseq-lrt-step-2"
  )
  parser$add_argument(
    "--threads",
    help = "Number of threads",
    type = "integer",
    default = 1
  )
  parser$add_argument(
    "--test_mode",
    help = "Run in test mode (first 500 rows only)",
    action = "store_true",
    default = FALSE
  )
  
  # Parse arguments
  args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
  
  # Convert boolean string values if needed
  for (arg_name in c("use_lfc_thresh", "test_mode")) {
    if (!is.null(args[[arg_name]])) {
      args[[arg_name]] <- convert_to_boolean(args[[arg_name]], FALSE)
    }
  }
  
  return(args)
} 