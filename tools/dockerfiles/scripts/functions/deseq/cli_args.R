#!/usr/bin/env Rscript
#
# Command-line argument handling for DESeq/DESeq2 differential expression analysis
#
# This file contains functions for parsing and validating command-line arguments
# for the main DESeq analysis workflow.
#
# Version: 0.1.0

# Load required libraries
suppressMessages(library(argparse))

#' Assert and validate command line arguments
#'
#' @param args The parsed arguments from ArgumentParser
#' @return Modified args with validated and processed values
assert_args <- function(args) {
  log_message("Checking input parameters")
  
  # Process aliases if not provided
  if (is.null(args$ualias) | is.null(args$talias)) {
    log_message("--ualias or --talias were not set, using default values based on expression file names")
    
    args$ualias <- character(0)
    for (i in 1:length(args$untreated)) {
      args$ualias <- append(args$ualias, head(unlist(
        strsplit(basename(args$untreated[i]), ".", fixed = TRUE)
      ), 1))
    }
    
    args$talias <- character(0)
    for (i in 1:length(args$treated)) {
      args$talias <- append(args$talias, head(unlist(
        strsplit(basename(args$treated[i]), ".", fixed = TRUE)
      ), 1))
    }
  } else {
    # Verify correct number of aliases
    if ((length(args$ualias) != length(args$untreated)) |
        (length(args$talias) != length(args$treated))) {
      log_error("Not correct number of inputs provided as -u, -t, -ua, -ut")
      quit(save = "no", status = 1, runLast = FALSE)
    }
  }

  # Check for minimum file requirements
  if (length(args$treated) == 1 || length(args$untreated) == 1) {
    log_warning("Only one file in a group. DESeq2 requires at least two replicates for accurate analysis.")
    args$batchfile <- NULL # reset batchfile to NULL. We don't need it for DESeq even if it was provided
  }

  # Process batch file if provided
  if (!is.null(args$batchfile)) {
    batch_metadata <- with_error_handling({
      read.table(
        args$batchfile,
        sep = get_file_type(args$batchfile),
        row.names = 1,
        col.names = c("name", "batch"),
        header = FALSE,
        stringsAsFactors = FALSE
      )
    })
    
    if (is.null(batch_metadata)) {
      log_error("Failed to read batch metadata file")
      args$batchfile <- NULL
      return(args)
    }
    
    log_message("Loaded batch metadata")
    rownames(batch_metadata) <- gsub("'|\"| ", "_", rownames(batch_metadata))
    
    if (all(is.element(c(args$ualias, args$talias), rownames(batch_metadata)))) {
      args$batchfile <- batch_metadata # dataframe
    } else {
      log_warning("Missing values in batch metadata file. Skipping multi-factor analysis")
      log_debug(paste("Expected:", paste(c(args$ualias, args$talias), collapse=", ")))
      log_debug(paste("Found:", paste(rownames(batch_metadata), collapse=", ")))
      args$batchfile <- NULL
    }
  }

  return(args)
}

#' Parse command line arguments for DESeq analysis
#'
#' @return Parsed and validated argument list
get_args <- function() {
  parser <- ArgumentParser(description = "Run DESeq/DESeq2 for untreated-vs-treated groups (condition-1-vs-condition-2)")
  
  # Input file parameters
  parser$add_argument(
    "-u", "--untreated",
    help = "Untreated (condition 1) CSV/TSV isoforms expression files",
    type = "character",
    required = TRUE,
    nargs = "+"
  )
  parser$add_argument(
    "-t", "--treated",
    help = "Treated (condition 2) CSV/TSV isoforms expression files",
    type = "character",
    required = TRUE,
    nargs = "+"
  )
  parser$add_argument(
    "-ua", "--ualias",
    help = "Unique aliases for untreated (condition 1) expression files. Default: basenames of -u without extensions",
    type = "character",
    nargs = "*"
  )
  parser$add_argument(
    "-ta", "--talias",
    help = "Unique aliases for treated (condition 2) expression files. Default: basenames of -t without extensions",
    type = "character",
    nargs = "*"
  )
  
  # Condition naming parameters
  parser$add_argument(
    "-un", "--uname",
    help = "Name for untreated (condition 1), use only letters and numbers",
    type = "character",
    default = "untreated"
  )
  parser$add_argument(
    "-tn", "--tname",
    help = "Name for treated (condition 2), use only letters and numbers",
    type = "character",
    default = "treated"
  )
  
  # Batch correction parameters
  parser$add_argument(
    "-bf", "--batchfile",
    help = paste(
      "Metadata file for multi-factor analysis. Headerless TSV/CSV file.",
      "First column - names from --ualias and --talias, second column - batch group name.",
      "Default: None"
    ),
    type = "character"
  )
  parser$add_argument(
    "--batchcorrection",
    help = paste(
      "Specifies the batch correction method to be applied.",
      "- 'combatseq' applies ComBat_seq at the beginning of the analysis, removing batch effects from the design formula before differential expression analysis.",
      "- 'model' applies removeBatchEffect from the limma package after differential expression analysis, incorporating batch effects into the model during DE analysis.",
      "- Default: none"
    ),
    type = "character",
    choices = c("none", "combatseq", "model"),
    default = "none"
  )
  
  # Statistical and filtering parameters
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
    "--rpkm_cutoff",
    help = paste(
      "RPKM cutoff for filtering genes. Genes with RPKM values below this threshold will be excluded from the analysis.",
      "Default: NULL (no filtering)"
    ),
    type = "integer",
    default = NULL
  )
  parser$add_argument(
    "--regulation",
    help = paste(
      "Direction of differential expression comparison. β is the log2 fold change.",
      "'both' for both up and downregulated genes (|β| > lfcThreshold for greaterAbs and |β| < lfcThreshold for lessAbs, with p-values being two-tailed or maximum of the upper and lower tests, respectively); ",
      "'up' for upregulated genes (β > lfcThreshold in condition2 compared to condition1); ",
      "'down' for downregulated genes (β < -lfcThreshold in condition2 compared to condition1). ",
      "Default: both"
    ),
    type = "character",
    choices = c("both", "up", "down"),
    default = "both"
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
  
  # Clustering parameters
  parser$add_argument(
    "--cluster",
    help = paste(
      "Hopach clustering method to be run on normalized read counts for the",
      "exploratory visualization part of the analysis. Default: do not run",
      "clustering"
    ),
    type = "character",
    choices = c("row", "column", "both")
  )
  parser$add_argument(
    "--scaling_type",
    help = paste(
      "Specifies the type of scaling to be applied to the expression data.",
      "- 'minmax' applies Min-Max scaling, normalizing values to a range of [-2, 2].",
      "- 'zscore' applies Z-score standardization, centering data to mean = 0 and standard deviation = 1.",
      "- Default: zscore"
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
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor")
  )
  parser$add_argument(
    "--columndist",
    help = paste(
      "Distance metric for HOPACH column clustering. Ignored if --cluster is not",
      "provided. Default: euclid"
    ),
    type = "character",
    default = "euclid",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor")
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
  
  # Output parameters
  parser$add_argument(
    "-o", "--output",
    help = "Output prefix. Default: deseq",
    type = "character",
    default = "./deseq"
  )
  parser$add_argument(
    "-d", "--digits",
    help = "Precision, number of digits to print. Default: 3",
    type = "integer",
    default = 3
  )
  parser$add_argument(
    "-p", "--threads",
    help = "Number of threads to use for parallel processing. Default: 1",
    type = "integer",
    default = 1
  )
  
  # Parse and clean arguments
  parsed_args <- parser$parse_args(gsub("'|\"| ", "_", commandArgs(trailingOnly = TRUE)))
  
  # Validate and process arguments
  validated_args <- assert_args(parsed_args)
  
  return(validated_args)
} 