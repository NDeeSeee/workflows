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
    help = "Use lfcthreshold as the null hypothesis value in the results function call"
  )
  
  parser$add_argument(
    "--rpkm_cutoff",
    type = "integer",
    default = NULL,
    help = "Integer cutoff for filtering rows in the expression data"
  )
  
  # Using the names directly as in CWL
  parser$add_argument(
    "--cluster",
    default = "none",
    choices = c("row", "column", "both", "none"),
    help = "Hopach clustering method to be run on normalized read counts"
  )
  
  parser$add_argument(
    "--rowdist",
    default = "cosangle",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor"),
    help = "Distance metric for HOPACH row clustering"
  )
  
  parser$add_argument(
    "--columndist",
    default = "euclid",
    choices = c("cosangle", "abscosangle", "euclid", "cor", "abscor"),
    help = "Distance metric for HOPACH column clustering"
  )
  
  parser$add_argument(
    "--k",
    type = "integer",
    default = 3,
    help = "Number of levels (depth) for Hopach clustering: min - 1, max - 15"
  )
  
  parser$add_argument(
    "--kmax",
    type = "integer",
    default = 5,
    help = "Maximum number of clusters at each level for Hopach clustering: min - 2, max - 9"
  )
  
  # Output arguments
  parser$add_argument(
    "--output",
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
          # Instead of filtering based on a specific pattern, accept all potential sample names
          # since they come from CWL and should be trusted
          args$name <- c(args$name, name_args)
          message(paste("Adding", length(name_args), "sample names from positional args"))
        }
      }
      
      # Additional check for name-input mismatch
      if (is.character(args$name) && length(args$name) == 1 && length(args$input) > 1) {
        # Special handling: if we have multiple input files but only one name value,
        # check if the name value could actually be a comma-separated list of names
        potential_names <- unlist(strsplit(args$name, ","))
        if (length(potential_names) > 1) {
          message("Detected comma-separated list of names, expanding to array")
          args$name <- potential_names
        } else {
          # If there's just one name and many inputs, but the single name isn't a comma-separated list,
          # check if there are any other arguments that could be sample names
          message("WARNING: Only one sample name found for multiple input files")
          message("Looking for additional sample names in remaining arguments...")
          
          # Try to detect sample names from the remaining arguments
          # This is more permissive than before
          remaining_args <- setdiff(all_args, c(flag_indices, flag_values, args$input))
          if (length(remaining_args) > 0) {
            potential_names <- remaining_args[!grepl("^--", remaining_args)]
            if (length(potential_names) > 0) {
              message(paste("Found", length(potential_names), "potential additional sample names"))
              args$name <- c(args$name, potential_names)
            }
          }
        }
      }
      
      # Final verification: ensure we have the same number of names as inputs
      if (length(args$input) > 0 && length(args$name) > 0) {
        if (length(args$name) == 1 && length(args$input) > 1) {
          # If we still have a mismatch, try to generate names from the input file paths
          message("WARNING: Input/name count mismatch. Generating sample names from input files.")
          args$name <- basename(args$input)
          args$name <- gsub("\\.(tsv|csv)$", "", args$name)
        } else if (length(args$name) < length(args$input)) {
          message("WARNING: Fewer sample names than input files. Some names may be missing.")
        }
      }
      
      # Parse numeric values
      for (arg_name in c("fdr", "lfcthreshold", "rpkm_cutoff", "k", "kmax", "threads")) {
        if (!is.null(args[[arg_name]]) && is.character(args[[arg_name]])) {
          if (grepl("^[0-9.]+$", args[[arg_name]])) {
            args[[arg_name]] <- as.numeric(args[[arg_name]])
          }
        }
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
  
  # Coerce args to a list if necessary (after tryCatch)
  if (!is.list(args)) {
    args <- as.list(args)
  }
  
  # Validate arguments
  args <- assert_args(args)
  
  # Coerce args to a list if necessary (after assert_args)
  if (!is.list(args)) {
    args <- as.list(args)
  }
  
  # Trim whitespace from name values
  if (!is.null(args$name) && length(args$name) > 0) {
    args$name <- trimws(args$name)
    message("Trimmed whitespace from sample names")
  }
  
  # Convert boolean string values to actual booleans if they came as strings
  for (arg_name in c("use_lfc_thresh", "lrt_only_mode", "test_mode")) {
    if (!is.null(args[[arg_name]])) {
      args[[arg_name]] <- convert_to_boolean(args[[arg_name]], FALSE)
    }
  }
  
  return(args)
} 