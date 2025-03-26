#!/usr/bin/env Rscript

# --- Command line argument parsing functions ---

# Function to parse command line arguments
get_args <- function() {
  # Get raw command line args for backup
  raw_args <- commandArgs(trailingOnly = TRUE)
  
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
  
  # Parse arguments safely with error handling
  tryCatch({
    args <- parser$parse_args()
  }, error = function(e) {
    # Check for unrecognized arguments error
    if (grepl("unrecognized arguments:", e$message) || grepl("expected one argument", e$message)) {
      message("Warning: Argument parsing error. Attempting to handle arguments manually.")
      
      # Get all command line arguments
      all_args <- commandArgs(trailingOnly = TRUE)
      
      # Initialize an empty list for our parsed arguments
      args <- list()
      
      # Initialize arrays for multi-value arguments
      array_args <- c("input", "name")
      for (arg in array_args) {
        args[[arg]] <- character(0)
      }
      
      # Make sure to explicitly extract required arguments first
      required_args <- c("meta", "design", "reduced")
      for (req_arg in required_args) {
        req_flag <- paste0("--", req_arg)
        arg_idx <- which(all_args == req_flag)
        if (length(arg_idx) > 0 && arg_idx[1] < length(all_args)) {
          args[[req_arg]] <- all_args[arg_idx[1] + 1]
          message(paste("Directly extracted required argument:", req_arg, "=", args[[req_arg]]))
        }
      }
      
      # First pass: process all flags with values
      i <- 1
      while (i <= length(all_args)) {
        current_arg <- all_args[i]
        
        # Check if this is a flag argument (starts with --)
        if (grepl("^--", current_arg)) {
          arg_name <- sub("^--", "", current_arg)
          
          # Check if the next item exists and is not a flag
          if (i < length(all_args) && !grepl("^--", all_args[i + 1])) {
            arg_value <- all_args[i + 1]
            
            # Handle array arguments - we need to append values
            if (arg_name %in% array_args) {
              args[[arg_name]] <- c(args[[arg_name]], arg_value)
            } else {
              # For scalar arguments, just set the value
              args[[arg_name]] <- arg_value
            }
            
            i <- i + 2  # Skip the value
          } else {
            # This is a boolean flag
            args[[arg_name]] <- TRUE
            i <- i + 1
          }
        } else {
          # This is a positional argument, classify it later
          i <- i + 1
        }
      }
      
      # At this point we have extracted all named parameters
      # Now collect all the --input and --name values specifically
      input_flags <- which(all_args == "--input")
      name_flags <- which(all_args == "--name")
      
      # Clear any incorrectly assigned values that might have been added to args$name
      if ("name" %in% names(args)) {
        args$name <- character(0)
      }
      
      # Extract input values correctly
      if (length(input_flags) > 0) {
        args$input <- character(0)
        for (idx in input_flags) {
          if (idx < length(all_args) && !grepl("^--", all_args[idx + 1])) {
            args$input <- c(args$input, all_args[idx + 1])
          }
        }
      }
      
      # Extract name values correctly
      if (length(name_flags) > 0) {
        args$name <- character(0)
        for (idx in name_flags) {
          if (idx < length(all_args) && !grepl("^--", all_args[idx + 1])) {
            args$name <- c(args$name, all_args[idx + 1])
          }
        }
      }
      
      # Second pass: Find any positional arguments
      positional_args <- all_args[!grepl("^--", all_args) & !grepl("^-[a-zA-Z]", all_args)]
      
      # Remove positional args that are values of flags
      flag_indices <- which(grepl("^--", all_args) | grepl("^-[a-zA-Z]", all_args))
      value_indices <- flag_indices + 1
      value_indices <- value_indices[value_indices <= length(all_args)]
      flag_values <- all_args[value_indices]
      positional_args <- setdiff(positional_args, flag_values)
      
      # If no positional args, return what we have
      if (length(positional_args) == 0) {
        message("No positional arguments found")
        return(args)
      }
      
      # Add to array arguments if we have any positional args
      if (length(positional_args) > 0) {
        message(paste("Found", length(positional_args), "potential positional arguments"))
        
        # Detect file paths and sample names
        file_pattern <- "\\.(tsv|csv)$"
        file_args <- positional_args[grepl(file_pattern, positional_args)]
        name_args <- positional_args[!grepl(file_pattern, positional_args)]
        
        if (length(file_args) > 0) {
          args$input <- c(args$input, file_args)
          message(paste("Added", length(file_args), "file paths to --input"))
        }
        
        # Only add names that match the sample naming pattern
        if (length(name_args) > 0) {
          # Filter only actual sample names (assume they start with ABC)
          sample_pattern <- "^ABSK\\d+"
          valid_names <- name_args[grepl(sample_pattern, name_args)]
          
          if (length(valid_names) > 0) {
            args$name <- c(args$name, valid_names)
            message(paste("Adding", length(valid_names), "sample names from positional args"))
          }
        }
      }
      
      # Parse numeric values
      for (arg_name in c("fdr", "lfcthreshold", "rpkm_cutoff", "k_hopach", "kmax_hopach", "threads")) {
        if (!is.null(args[[arg_name]]) && is.character(args[[arg_name]])) {
          if (grepl("^[0-9.]+$", args[[arg_name]])) {
            args[[arg_name]] <- as.numeric(args[[arg_name]])
          }
        }
      }
      
      # Convert boolean string values
      for (arg_name in c("use_lfc_thresh", "lrt_only_mode", "test_mode")) {
        if (!is.null(args[[arg_name]]) && is.character(args[[arg_name]])) {
          if (tolower(args[[arg_name]]) %in% c("true", "t", "yes", "y", "1")) {
            args[[arg_name]] <- TRUE
          } else if (tolower(args[[arg_name]]) %in% c("false", "f", "no", "n", "0")) {
            args[[arg_name]] <- FALSE
          }
        }
      }
      
      # Set up parameter aliases for compatibility
      if (!is.null(args$k)) {
        args$k_hopach <- args$k
      }
      if (!is.null(args$kmax)) {
        args$kmax_hopach <- args$kmax
      }
      if (!is.null(args$cluster)) {
        args$cluster_method <- args$cluster
      }
      if (!is.null(args$rowdist)) {
        args$row_distance <- args$rowdist
      }
      if (!is.null(args$columndist)) {
        args$column_distance <- args$columndist
      }
      if (!is.null(args$output)) {
        args$output_prefix <- args$output
      }
      
      # Show what we parsed
      message("Manually parsed arguments:")
      for (arg_name in names(args)) {
        if (length(args[[arg_name]]) > 1) {
          message(paste0("  ", arg_name, ": [", paste(head(args[[arg_name]], 3), collapse=", "), 
                        if(length(args[[arg_name]]) > 3) "..." else "", "] (", length(args[[arg_name]]), " items)"))
        } else {
          message(paste0("  ", arg_name, ": ", args[[arg_name]]))
        }
      }
      
      return(args)
    } else {
      # For other errors, just stop with the error message
      stop(e$message)
    }
  })
  
  # Validate arguments
  args <- assert_args(args)
  
  # Trim whitespace from name values
  if (!is.null(args$name) && length(args$name) > 0) {
    args$name <- trimws(args$name)
    message("Trimmed whitespace from sample names")
  }
  
  # Final check for required arguments before returning
  required_args <- c("meta", "design", "reduced")
  for (req_arg in required_args) {
    if (is.null(args[[req_arg]]) || !req_arg %in% names(args)) {
      # One last attempt to recover from raw command line arguments
      raw_args <- commandArgs(trailingOnly = TRUE)
      arg_flag <- paste0("--", req_arg)
      arg_idx <- which(raw_args == arg_flag)
      if (length(arg_idx) > 0 && arg_idx[1] < length(raw_args)) {
        args[[req_arg]] <- raw_args[arg_idx[1] + 1]
        message(paste("Final recovery of required argument:", req_arg, "=", args[[req_arg]]))
      } else {
        message(paste("Error: Could not recover required argument:", req_arg))
      }
    }
  }
  
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
  # Check for required arguments
  required_args <- c("meta", "design", "reduced")
  missing_required <- required_args[!required_args %in% names(args)]
  
  if (length(missing_required) > 0) {
    # Before failing, check if any required arguments may have been parsed 
    # correctly but are not accessible due to variable scope issues
    # Try to extract these from the original command arguments if possible
    all_args <- commandArgs(trailingOnly = TRUE)
    
    for (missing_arg in missing_required) {
      arg_flag <- paste0("--", missing_arg)
      arg_idx <- which(all_args == arg_flag)
      if (length(arg_idx) > 0 && arg_idx[1] < length(all_args)) {
        args[[missing_arg]] <- all_args[arg_idx[1] + 1]
        message(paste("Recovered missing required argument:", missing_arg, "=", args[[missing_arg]]))
        # Remove from missing list
        missing_required <- setdiff(missing_required, missing_arg)
      }
    }
    
    # If we still have missing arguments, abort
    if (length(missing_required) > 0) {
      message(paste("Error: Missing required arguments:", paste(missing_required, collapse=", ")))
      quit(save = "no", status = 1, runLast = FALSE)
    }
  }
  
  # Create empty arrays for input/name if they don't exist
  if (is.null(args$input)) args$input <- character(0)
  if (is.null(args$name)) args$name <- character(0)
  
  # Validate input and name parameters
  if (length(args$input) == 0 || length(args$name) == 0) {
    message("Error: --input and --name parameters are required")
    message("Input files found: ", length(args$input))
    message("Sample names found: ", length(args$name))
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  if (length(args$input) != length(args$name)) {
    message("Error: --input and --name have different number of values")
    message(paste("Number of input files:", length(args$input)))
    message(paste("Number of sample names:", length(args$name)))
    message("Input files:")
    print(args$input)
    message("Sample names:")
    print(args$name)
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  # Validate design formula
  tryCatch(
    expr = {
      # Try to load design formula
      design_formula <- as.formula(args$design)
    },
    error = function(e) {
      message(paste0("Error: failed to load --design ", args$design, " as formula"))
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
      message(paste0("Error: failed to load --reduced ", args$reduced, " as formula"))
      quit(save = "no", status = 1, runLast = FALSE)
    }
  )
  
  # Check if metadata file exists
  if (!file.exists(args$meta)) {
    message(paste("Error: Metadata file does not exist:", args$meta))
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  # Validate batch correction
  if (!is.null(args$batchcorrection) && args$batchcorrection != "none") {
    # Try to read the metadata file
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
        message(paste0("Error: failed to read metadata file: ", e$message))
        quit(save = "no", status = 1, runLast = FALSE)
      }
    )
    
    if (!"batch" %in% colnames(meta_data)) {
      message("Warning: batch correction requested but 'batch' column not found in metadata file")
      message("Proceeding without batch correction")
      args$batchcorrection <- "none"
    } else if (!is.numeric(meta_data$batch)) {
      message("Warning: 'batch' column in metadata file is not numeric")
      message("Converting batch to numeric values")
      # No need to quit here, we'll handle conversion elsewhere
    }
  }
  
  # Validate rpkm_cutoff (if provided)
  if (!is.null(args$rpkm_cutoff)) {
    if (!is.numeric(args$rpkm_cutoff) || args$rpkm_cutoff < 0) {
      message("Warning: --rpkm_cutoff must be a non-negative integer, using default NULL")
      args$rpkm_cutoff <- NULL
    }
  }
  
  # Validate k_hopach parameter (1-15)
  if (!is.null(args$k_hopach) && (args$k_hopach < 1 || args$k_hopach > 15)) {
    message("Warning: --k_hopach must be between 1 and 15, using default value 3")
    args$k_hopach <- 3
  }
  
  # Validate kmax_hopach parameter (2-9)
  if (!is.null(args$kmax_hopach) && (args$kmax_hopach < 2 || args$kmax_hopach > 9)) {
    message("Warning: --kmax_hopach must be between 2 and 9, using default value 5")
    args$kmax_hopach <- 5
  }
  
  # Validate fdr parameter
  if (!is.null(args$fdr) && (args$fdr < 0 || args$fdr > 1)) {
    message("Warning: --fdr must be between 0 and 1, using default value 0.1")
    args$fdr <- 0.1
  }
  
  # Set default values for any missing parameters
  if (is.null(args$batchcorrection)) args$batchcorrection <- "none"
  if (is.null(args$fdr)) args$fdr <- 0.1
  if (is.null(args$lfcthreshold)) args$lfcthreshold <- 0.59
  if (is.null(args$use_lfc_thresh)) args$use_lfc_thresh <- FALSE
  if (is.null(args$cluster_method)) args$cluster_method <- "none"
  if (is.null(args$row_distance)) args$row_distance <- "cosangle"
  if (is.null(args$column_distance)) args$column_distance <- "euclid"
  if (is.null(args$k_hopach)) args$k_hopach <- 3
  if (is.null(args$kmax_hopach)) args$kmax_hopach <- 5
  if (is.null(args$output_prefix)) args$output_prefix <- "./deseq_lrt_step_1"
  if (is.null(args$threads)) args$threads <- 1
  if (is.null(args$lrt_only_mode)) args$lrt_only_mode <- FALSE
  if (is.null(args$test_mode)) args$test_mode <- FALSE
  
  return(args)
}

# This namespace contains exported parameter handling functions
params <- new.env()

#' Get command line arguments with standard parsing
#' @export
params$get_cli_args <- function() {
  return(get_args())
}

#' Assert that required arguments are present and valid
#' @param args List of arguments to validate
#' @export
params$assert_args <- function(args) {
  # Check required arguments are present
  required_args <- c("meta", "design", "reduced", "input")
  missing_args <- required_args[!required_args %in% names(args)]
  
  if (length(missing_args) > 0) {
    # Before failing, try to recover missing arguments from command line
    all_args <- commandArgs(trailingOnly = TRUE)
    
    for (missing_arg in missing_args) {
      arg_flag <- paste0("--", missing_arg)
      arg_idx <- which(all_args == arg_flag)
      if (length(arg_idx) > 0 && arg_idx[1] < length(all_args)) {
        args[[missing_arg]] <- all_args[arg_idx[1] + 1]
        message(paste("Recovered missing required argument:", missing_arg, "=", args[[missing_arg]]))
        # Remove from missing list
        missing_args <- setdiff(missing_args, missing_arg)
      }
    }
    
    # If we still have missing arguments, stop
    if (length(missing_args) > 0) {
      stop(paste("Missing required arguments:", paste(missing_args, collapse=", ")))
    }
  }
  
  # Check that input files exist
  for (input_file in args$input) {
    if (!file.exists(input_file)) {
      stop(paste("Input file does not exist:", input_file))
    }
  }
  
  # Check that metadata file exists
  if (!file.exists(args$meta)) {
    stop(paste("Metadata file does not exist:", args$meta))
  }
  
  # Validate design formulas
  tryCatch({
    as.formula(args$design)
    as.formula(args$reduced)
  }, error = function(e) {
    stop(paste("Invalid formula:", e$message))
  })
  
  return(invisible(TRUE))
}

# Export the params namespace
assign("params", params, envir = .GlobalEnv) 