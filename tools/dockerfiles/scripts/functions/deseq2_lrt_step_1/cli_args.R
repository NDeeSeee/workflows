#!/usr/bin/env Rscript

# --- Command line argument parsing functions ---

# Function to parse command line arguments
get_args <- function() {
  parser <- argparse::ArgumentParser(
    description = "Run DESeq2 analysis with Likelihood Ratio Test (LRT)",
    formatter_class = "argparse.ArgumentDefaultsHelpFormatter"
  )
  
  # Input data arguments
  parser$add_argument(
    "--input",
    action = "append",
    help = "List of input files with expression data (CSV or TSV format)"
  )
  parser$add_argument(
    "--name",
    action = "append",
    help = "Names for input files (in the same order as input files)"
  )
  parser$add_argument(
    "--meta",
    required = TRUE,
    help = "Metadata file in CSV or TSV format"
  )
  
  # Analysis parameters
  parser$add_argument(
    "--design",
    required = TRUE,
    help = "Design formula for DESeq2 (e.g., '~condition+batch')"
  )
  parser$add_argument(
    "--reduced",
    required = TRUE,
    help = "Reduced design formula for LRT (e.g., '~batch')"
  )
  
  # CWL-aligned parameters
  parser$add_argument(
    "--batchcorrection",
    default = "none",
    choices = c("none", "combatseq", "model"),
    help = "Batch correction method: 'none', 'combatseq', or 'model'"
  )
  
  parser$add_argument(
    "--scaling_type",
    default = "zscore",
    choices = c("minmax", "zscore"),
    help = "Scaling type for expression data: 'minmax' or 'zscore'"
  )
  
  parser$add_argument(
    "--fdr",
    type = "double",
    default = 0.1,
    help = "FDR threshold for significance"
  )
  
  parser$add_argument(
    "--lfcthreshold",
    type = "double",
    default = 0.59,
    help = "Log2 fold change threshold for determining significant differential expression"
  )
  
  parser$add_argument(
    "--use_lfc_thresh",
    action = "store_true",
    default = FALSE,
    help = "How to apply LFC threshold: TRUE - use in hypothesis testing (as null), FALSE - apply as post-filtering. Note: For LRT tests, LFC threshold is not used in testing."
  )
  
  parser$add_argument(
    "--rpkm_cutoff",
    type = "integer",
    default = NULL,
    help = "Integer cutoff for filtering rows in the expression data"
  )
  
  # Keep original CWL parameter names but store in the new parameter names
  # This avoids the ambiguous option error
  parser$add_argument(
    "--cluster",
    dest = "cluster_method",
    default = "none",
    choices = c("row", "column", "both", "none"),
    help = "Hopach clustering method to be run on normalized read counts"
  )
  
  parser$add_argument(
    "--rowdist",
    dest = "row_distance",
    default = "cosangle",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor"),
    help = "Distance metric for HOPACH row clustering"
  )
  
  parser$add_argument(
    "--columndist",
    dest = "column_distance",
    default = "euclid",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor"),
    help = "Distance metric for HOPACH column clustering"
  )
  
  parser$add_argument(
    "--k",
    dest = "k_hopach",
    type = "integer",
    default = 3,
    help = "Number of levels (depth) for Hopach clustering: min - 1, max - 15"
  )
  
  parser$add_argument(
    "--kmax",
    dest = "kmax_hopach",
    type = "integer",
    default = 5,
    help = "Maximum number of clusters at each level for Hopach clustering: min - 2, max - 9"
  )
  
  # Output arguments
  parser$add_argument(
    "--output",
    dest = "output_prefix",
    default = "./deseq_lrt_step_1",
    help = "Output prefix for generated files"
  )
  
  parser$add_argument(
    "--threads",
    type = "integer",
    default = 1,
    help = "Number of threads to use for parallel processing"
  )
  
  parser$add_argument(
    "--lrt_only_mode",
    action = "store_true",
    default = FALSE,
    help = "Run LRT only, no contrasts"
  )
  
  parser$add_argument(
    "--test_mode",
    action = "store_true",
    default = FALSE,
    help = "Run for test, only first 500 rows"
  )
  
  # Parse arguments
  args <- parser$parse_args()
  
  # Validate arguments
  args <- assert_args(args)
  
  return(args)
}

# Function to handle parameter naming compatibility
handle_parameter_compatibility <- function(args) {
  # NOTE: This function is no longer needed as we're handling parameter compatibility
  # directly in the ArgumentParser using the 'dest' parameter
  return(args)
}

# Function to validate command line arguments
assert_args <- function(args) {
  # Validate input and name parameters
  if (is.null(args$input) || is.null(args$name)) {
    print("Exiting: --input and --name parameters are required")
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  if (length(args$input) != length(args$name)) {
    print("Exiting: --input and --name have different number of values")
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  # Validate design formula
  tryCatch(
    expr = {
      # Try to load design formula
      design_formula <- as.formula(args$design)
    },
    error = function(e) {
      print(paste0("Exiting: failed to load --design ", args$design, " as formula"))
      quit(save = "no", status = 1, runLast = FALSE)
    }
  )
  
  # Validate reduced formula
  tryCatch(
    expr = {
      # Try to load reduced formula
      reduced_formula <- as.formula(args$reduced)
    },
    error = function(e) {
      print(paste0("Exiting: failed to load --reduced ", args$reduced, " as formula"))
      quit(save = "no", status = 1, runLast = FALSE)
    }
  )
  
  # Validate batch correction
  if (args$batchcorrection != "none") {
    # Check if metadata file exists
    if (!file.exists(args$meta)) {
      print("Exiting: metadata file does not exist")
      quit(save = "no", status = 1, runLast = FALSE)
    }
    
    # Check if 'batch' column exists in metadata
    meta_data <- tryCatch(
      {
        # Try to read the file
        if (endsWith(args$meta, ".csv")) {
          read.csv(args$meta, check.names = FALSE, stringsAsFactors = FALSE)
        } else if (endsWith(args$meta, ".tsv")) {
          read.delim(args$meta, check.names = FALSE, stringsAsFactors = FALSE)
        } else {
          # Default to CSV
          read.csv(args$meta, check.names = FALSE, stringsAsFactors = FALSE)
        }
      },
      error = function(e) {
        print(paste0("Exiting: failed to read metadata file: ", e$message))
        quit(save = "no", status = 1, runLast = FALSE)
      }
    )
    
    if (!"batch" %in% colnames(meta_data)) {
      print("Exiting: batch correction requested but 'batch' column not found in metadata file")
      quit(save = "no", status = 1, runLast = FALSE)
    }
    
    # Check if batch column is numeric
    if (!is.numeric(meta_data$batch)) {
      print("Exiting: 'batch' column in metadata file must be numeric")
      quit(save = "no", status = 1, runLast = FALSE)
    }
  }
  
  # Validate rpkm_cutoff (if provided)
  if (!is.null(args$rpkm_cutoff)) {
    if (!is.numeric(args$rpkm_cutoff) || args$rpkm_cutoff < 0) {
      print("Exiting: --rpkm_cutoff must be a non-negative integer")
      quit(save = "no", status = 1, runLast = FALSE)
    }
  }
  
  # Validate k_hopach parameter (1-15)
  if (args$k_hopach < 1 || args$k_hopach > 15) {
    print("Exiting: --k_hopach must be between 1 and 15")
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  # Validate kmax_hopach parameter (2-9)
  if (args$kmax_hopach < 2 || args$kmax_hopach > 9) {
    print("Exiting: --kmax_hopach must be between 2 and 9")
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  # Validate fdr parameter
  if (args$fdr < 0 || args$fdr > 1) {
    print("Exiting: --fdr must be between 0 and 1")
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  return(args)
} 