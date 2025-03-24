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
  parser$add_argument(
    "--rpkm-cutoff",
    type = "double",
    default = NULL,
    help = "RPKM cutoff for filtering genes (NULL for no filtering)"
  )
  parser$add_argument(
    "--mincounts",
    type = "integer",
    default = 10,
    help = "Minimum count across all samples for filtering genes"
  )
  parser$add_argument(
    "--fdr",
    type = "double",
    default = 0.1,
    help = "FDR threshold for significance"
  )
  parser$add_argument(
    "--batchcorrection",
    default = "none",
    help = "Batch correction method: 'none', 'limma', 'combat', or 'combatseq'"
  )
  
  # Output arguments
  parser$add_argument(
    "--output",
    default = "output",
    help = "Output directory path"
  )
  parser$add_argument(
    "--prefix",
    default = "deseq2_lrt",
    help = "Prefix for output files"
  )
  
  # Optional debugging arguments
  parser$add_argument(
    "--debug",
    action = "store_true",
    default = FALSE,
    help = "Enable debug mode with additional logging"
  )
  parser$add_argument(
    "--verbose",
    action = "store_true",
    default = FALSE,
    help = "Enable verbose error messages"
  )
  parser$add_argument(
    "--clean-logs",
    action = "store_true",
    default = TRUE,
    help = "Clean previous log file on start"
  )
  parser$add_argument(
    "--threads",
    type = "integer",
    default = 4,
    help = "Number of threads to use for parallel processing"
  )
  parser$add_argument(
    "--test-mode",
    action = "store_true",
    default = FALSE,
    help = "Run in test mode with limited data (for development only)"
  )
  
  # Parse arguments
  args <- parser$parse_args()
  
  # Validate arguments
  args <- assert_args(args)
  
  return(args)
}

# Function to validate command line arguments
assert_args <- function(args) {
  if (length(args$input) != length(args$name)) {
    print("Exiting: --input and --name have different number of values")
    quit(
      save = "no",
      status = 1,
      runLast = FALSE
    )
  }
  tryCatch(
    expr = {
      # Try to load design formula
      design_formula <- as.formula(tolower(args$design))
    },
    error = function(e) {
      print(paste0("Exiting: failed to load --design ", tolower(args$design), " as formula"))
      quit(
        save = "no",
        status = 1,
        runLast = FALSE
      )
    }
  )
  tryCatch(
    expr = {
      # Try to load reduced formula
      reduced_formula <- as.formula(tolower(args$reduced))
    },
    error = function(e) {
      print(paste0("Exiting: failed to load --reduced ", tolower(args$reduced), " as formula"))
      quit(
        save = "no",
        status = 1,
        runLast = FALSE
      )
    }
  )
  return(args)
} 