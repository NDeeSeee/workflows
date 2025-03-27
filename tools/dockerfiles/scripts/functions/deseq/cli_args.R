#!/usr/bin/env Rscript
#
# Command-line argument handling for DESeq/DESeq2 differential expression analysis
#
# This file contains functions for parsing and validating command-line arguments
# for the main DESeq analysis workflow.
#
# Version: 0.1.0

#' Assert and validate command line arguments
#'
#' @param args The parsed arguments from ArgumentParser
#' @return Modified args with validated and processed values
assert_args <- function(args) {
  log_message("Checking input parameters")
  
  # Process aliases if not provided
  if (is.null(args$untreated_sample_names) | is.null(args$treated_sample_names)) {
    log_message("--untreated_sample_names or --treated_sample_names were not set, using default values based on expression file names")
    
    args$untreated_sample_names <- character(0)
    for (i in 1:length(args$untreated_files)) {
      args$untreated_sample_names <- append(args$untreated_sample_names, head(unlist(
        strsplit(basename(args$untreated_files[i]), ".", fixed = TRUE)
      ), 1))
    }
    
    args$treated_sample_names <- character(0)
    for (i in 1:length(args$treated_files)) {
      args$treated_sample_names <- append(args$treated_sample_names, head(unlist(
        strsplit(basename(args$treated_files[i]), ".", fixed = TRUE)
      ), 1))
    }
  } else {
    # Verify correct number of aliases
    if ((length(args$untreated_sample_names) != length(args$untreated_files)) |
        (length(args$treated_sample_names) != length(args$treated_files))) {
      log_error("Not correct number of inputs provided for files and sample names")
      quit(save = "no", status = 1, runLast = FALSE)
    }
  }

  # Check for minimum file requirements
  if (length(args$treated_files) == 1 || length(args$untreated_files) == 1) {
    log_warning("Only one file in a group. DESeq2 requires at least two replicates for accurate analysis.")
    args$batch_file <- NULL # reset batch_file to NULL. We don't need it for DESeq even if it was provided
  }

  # Process batch file if provided
  if (!is.null(args$batch_file)) {
    batch_metadata <- with_error_handling({
      read.table(
        args$batch_file,
        sep = get_file_type(args$batch_file),
        row.names = 1,
        col.names = c("name", "batch"),
        header = FALSE,
        stringsAsFactors = FALSE
      )
    })
    
    if (is.null(batch_metadata)) {
      log_error("Failed to read batch metadata file")
      args$batch_file <- NULL
      return(args)
    }
    
    log_message("Loaded batch metadata")
    rownames(batch_metadata) <- gsub("'|\"| ", "_", rownames(batch_metadata))
    
    if (all(is.element(c(args$untreated_sample_names, args$treated_sample_names), rownames(batch_metadata)))) {
      args$batch_file <- batch_metadata # dataframe
    } else {
      log_warning("Missing values in batch metadata file. Skipping multi-factor analysis")
      log_debug(paste("Expected:", paste(c(args$untreated_sample_names, args$treated_sample_names), collapse=", ")))
      log_debug(paste("Found:", paste(rownames(batch_metadata), collapse=", ")))
      args$batch_file <- NULL
    }
  }

  # Convert boolean string values if they came as strings
  for (arg_name in c("use_lfc_thresh")) {
    if (!is.null(args[[arg_name]])) {
      args[[arg_name]] <- convert_to_boolean(args[[arg_name]], FALSE)
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
    "-u", "--untreated_files",
    help = "Untreated (condition 1) CSV/TSV isoforms expression files",
    type = "character",
    required = TRUE,
    nargs = "+"
  )
  parser$add_argument(
    "-t", "--treated_files",
    help = "Treated (condition 2) CSV/TSV isoforms expression files",
    type = "character",
    required = TRUE,
    nargs = "+"
  )
  parser$add_argument(
    "-ua", "--untreated_sample_names",
    help = "Unique aliases for untreated (condition 1) expression files. Default: basenames of -u without extensions",
    type = "character",
    nargs = "*"
  )
  parser$add_argument(
    "-ta", "--treated_sample_names",
    help = "Unique aliases for treated (condition 2) expression files. Default: basenames of -t without extensions",
    type = "character",
    nargs = "*"
  )
  
  # Condition naming parameters
  parser$add_argument(
    "-un", "--untreated_name",
    help = "Name for untreated (condition 1), use only letters and numbers",
    type = "character",
    default = "untreated"
  )
  parser$add_argument(
    "-tn", "--treated_name",
    help = "Name for treated (condition 2), use only letters and numbers",
    type = "character",
    default = "treated"
  )
  
  # Batch correction parameters
  parser$add_argument(
    "-bf", "--batch_file",
    help = paste(
      "Metadata file for multi-factor analysis. Headerless TSV/CSV file.",
      "First column - names from --untreated_sample_names and --treated_sample_names, second column - batch group name.",
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
    action = "store_true",
    default = FALSE
  )
  
  # Clustering parameters
  parser$add_argument(
    "--cluster_method",
    help = paste(
      "Hopach clustering method to be run on normalized read counts for the",
      "exploratory visualization part of the analysis. Default: none"
    ),
    type = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
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
    "--row_distance",
    help = paste(
      "Distance metric for HOPACH row clustering. Ignored if --cluster_method is not",
      "provided. Default: cosangle"
    ),
    type = "character",
    default = "cosangle",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor")
  )
  parser$add_argument(
    "--column_distance",
    help = paste(
      "Distance metric for HOPACH column clustering. Ignored if --cluster_method is not",
      "provided. Default: euclid"
    ),
    type = "character",
    default = "euclid",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor")
  )
  parser$add_argument(
    "--k_hopach",
    help = "Number of levels (depth) for Hopach clustering: min - 1, max - 15. Default: 3.",
    type = "integer",
    default = 3
  )
  parser$add_argument(
    "--kmax_hopach",
    help = "Maximum number of clusters at each level for Hopach clustering: min - 2, max - 9. Default: 5.",
    type = "integer",
    default = 5
  )
  
  # Output parameters
  parser$add_argument(
    "-o", "--output_prefix",
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
  
  # Handle legacy parameter aliases
  if (is.null(parsed_args$untreated_sample_names) && !is.null(parsed_args$ualias)) {
    parsed_args$untreated_sample_names <- parsed_args$ualias
  }
  if (is.null(parsed_args$treated_sample_names) && !is.null(parsed_args$talias)) {
    parsed_args$treated_sample_names <- parsed_args$talias
  }
  if (is.null(parsed_args$untreated_files) && !is.null(parsed_args$untreated)) {
    parsed_args$untreated_files <- parsed_args$untreated
  }
  if (is.null(parsed_args$treated_files) && !is.null(parsed_args$treated)) {
    parsed_args$treated_files <- parsed_args$treated
  }
  if (is.null(parsed_args$batch_file) && !is.null(parsed_args$batchfile)) {
    parsed_args$batch_file <- parsed_args$batchfile
  }
  if (is.null(parsed_args$cluster_method) && !is.null(parsed_args$cluster)) {
    parsed_args$cluster_method <- parsed_args$cluster
  }
  if (is.null(parsed_args$row_distance) && !is.null(parsed_args$rowdist)) {
    parsed_args$row_distance <- parsed_args$rowdist
  }
  if (is.null(parsed_args$column_distance) && !is.null(parsed_args$columndist)) {
    parsed_args$column_distance <- parsed_args$columndist
  }
  if (is.null(parsed_args$output_prefix) && !is.null(parsed_args$output)) {
    parsed_args$output_prefix <- parsed_args$output
  }
  
  # Validate and process arguments
  validated_args <- assert_args(parsed_args)
  
  return(validated_args)
} 