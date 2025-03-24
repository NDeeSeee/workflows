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
    "-dsq",
    "--dsq_obj_data",
    help = "RDS file containing contrasts and expression data from step 1",
    type = "character",
    required = TRUE
  )
  parser$add_argument(
    "-cd",
    "--contrast_df",
    help = "TSV file containing contrasts data",
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
    "--rpkm_cutoff",
    help = paste(
      "RPKM cutoff for filtering genes. Genes with RPKM values below this threshold will be excluded from the analysis.",
      "Default: NULL (no filtering)"
    ),
    type = "integer",
    default = NULL
  )
  parser$add_argument(
    "--scaling_type",
    help = paste(
      "Specifies the type of scaling to be applied to the expression data.",
      "- 'minmax' applies Min-Max scaling, normalizing values to a range of [-2, 2].",
      "- 'zscore' applies Z-score standardization, centering data to mean = 0 and standard deviation = 1.",
      "- Default: none (no scaling applied)."
    ),
    type = "character",
    choices = c("minmax", "zscore"),
    default = "zscore"
  )
  parser$add_argument(
    "--rowdist",
    help = paste(
      "Distance metric for HOPACH row clustering. Ignored if --cluster is not",
      "provided. Default: cosangle"
    ),
    type = "character",
    default = "cosangle",
    choices = c(
      "cosangle",
      "euclid",
      "abseuclid",
      "cor",
      "abscor"
    )
  )
  parser$add_argument(
    "--columndist",
    help = paste(
      "Distance metric for HOPACH column clustering. Ignored if --cluster is not",
      "provided. Default: euclid"
    ),
    type = "character",
    default = "euclid",
    choices = c(
      "cosangle",
      "euclid",
      "abseuclid",
      "cor",
      "abscor"
    )
  )
  parser$add_argument(
    "--cluster",
    help = paste(
      "Hopach clustering method to be run on normalized read counts for the",
      "exploratory visualization part of the analysis. Default: do not run",
      "clustering"
    ),
    type = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
  )
  parser$add_argument(
    "--k",
    help = "Number of levels (depth) for Hopach clustering: min - 1, max - 15. Default: 3.",
    type = "integer",
    default = 3
  )
  parser$add_argument(
    "--kmax",
    help = "Maximum number of clusters at each level for Hopach clustering: min - 2, max - 9. Default: 5.",
    type = "integer",
    default = 5
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
    help = "Run for test, only first 500 rows, clustering minimised",
    action = "store_true",
    default = FALSE
  )
  args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
  return(args)
} 