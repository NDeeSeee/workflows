#!/usr/bin/env Rscript
#
# DESeq2 LRT Analysis - Step 1
#
# This script performs differential expression analysis using DESeq2 with 
# Likelihood Ratio Test (LRT). It's been refactored for better maintainability
# with functions organized into separate files.
#
# Version: 0.1.3

# Set options
options(warn = -1)
options(rlang_backtrace_on_error = "full")
options("width" = 400)
options(error = function() {
  message("An unexpected error occurred. Aborting script.")
  quit(save = "no", status = 1, runLast = FALSE)
})

# Set memory management options for large datasets
options(future.globals.maxSize = 4000 * 1024^2)  # 4GB max for global data
options(expressions = 5000)  # Increase expression stack size

# Configure garbage collection behavior
gcinfo(FALSE)  # Disable GC messages by default
options(gc.aggressiveness = 0)  # Default GC behavior

# Load required libraries
suppressMessages({
  # Core packages
  library(argparse)
  library(BiocParallel)
  library(DESeq2)
  
  # Data manipulation
  library(tidyverse)
  library(data.table)
  library(conflicted)  # For resolving namespace conflicts
  
  # Batch correction
  library(limma)
  
  # Visualization
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(ggrepel)
  library(plotly)
  
  # GCT export
  library(cmapR)
  
  # Utilities
  library(pryr)        # For memory usage tracking
  library(rlang)
  library(stringr)
  library(glue)
})

# Resolve namespace conflicts explicitly
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")
conflicted::conflict_prefer("slice", "dplyr")
conflicted::conflict_prefer("mutate", "dplyr")
conflicted::conflict_prefer("group_by", "dplyr")
conflicted::conflict_prefer("summarize", "dplyr")
conflicted::conflict_prefer("arrange", "dplyr")
conflicted::conflict_prefer("%>%", "magrittr")
conflicted::conflict_prefer("%in%", "base")

# Additional conflict resolutions
conflicted::conflict_prefer("intersect", "base")
conflicted::conflict_prefer("setdiff", "base")
conflicted::conflict_prefer("union", "base")
conflicted::conflict_prefer("as.data.frame", "base")
conflicted::conflict_prefer("lag", "stats")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("first", "dplyr")
conflicted::conflict_prefer("last", "dplyr")
conflicted::conflict_prefer("layout", "graphics")
conflicted::conflict_prefer("plot", "graphics")
conflicted::conflict_prefer("desc", "dplyr")

# Pre-process command line arguments to handle unexpected positional arguments
preprocess_args <- function() {
  # Get all arguments
  all_args <- commandArgs(trailingOnly = TRUE)
  
  # Print original arguments for debugging
  message("Original arguments:")
  for (i in seq_along(all_args)) {
    message(paste("  ", i, ":", all_args[i]))
  }
  
  # Check if we have any arguments at all
  if (length(all_args) == 0) {
    message("No command line arguments found!")
    return(all_args)
  }
  
  # For array-based arguments, we need to preserve them during reconstruction
  array_args <- c("--input", "--name")
  array_values <- list()
  
  # Collect all array arguments first
  for (arg in array_args) {
    indices <- which(all_args == arg)
    if (length(indices) > 0) {
      values <- c()
      for (idx in indices) {
        if (idx < length(all_args) && !grepl("^--", all_args[idx+1])) {
          values <- c(values, all_args[idx+1])
        }
      }
      if (length(values) > 0) {
        array_values[[arg]] <- values
      }
    }
  }
  
  # Extract key required arguments first (meta, design, reduced)
  required_args <- c("--meta", "--design", "--reduced")
  required_values <- list()
  for (req_arg in required_args) {
    arg_index <- which(all_args == req_arg)
    if (length(arg_index) > 0 && arg_index[1] < length(all_args)) {
      required_values[[req_arg]] <- all_args[arg_index[1] + 1]
    }
  }
  
  # Extract all other singleton args (non-array args)
  other_args <- list()
  i <- 1
  while (i <= length(all_args)) {
    if (grepl("^--", all_args[i]) && !(all_args[i] %in% array_args)) {
      # If next arg exists and is not a flag
      if (i < length(all_args) && !grepl("^--", all_args[i+1])) {
        other_args[[all_args[i]]] <- all_args[i+1]
        i <- i + 2
      } else {
        # Just a flag
        other_args[[all_args[i]]] <- TRUE
        i <- i + 1
      }
    } else {
      # Skip array args (already handled) and positional args (handled later)
      i <- i + 1
    }
  }
  
  # Check if we have any arguments without -- prefix (positional arguments)
  positional_args <- all_args[!grepl("^--", all_args) & !grepl("^-[a-zA-Z]", all_args)]
  
  # Remove positional args that are values of flags
  flag_indices <- which(grepl("^--", all_args) | grepl("^-[a-zA-Z]", all_args))
  value_indices <- flag_indices + 1
  value_indices <- value_indices[value_indices <= length(all_args)]
  positional_values <- all_args[value_indices]
  positional_args <- setdiff(positional_args, positional_values)
  
  # If we have positional arguments, try to classify them
  if (length(positional_args) > 0) {
    message("Found positional arguments:", paste(positional_args, collapse=", "))
    
    # Separate positional arguments into names and inputs based on file path pattern
    file_pattern <- "\\.(tsv|csv)$"
    input_files <- positional_args[grepl(file_pattern, positional_args)]
    sample_names <- positional_args[!grepl(file_pattern, positional_args)]
    
    # If we have files or names, add them to our array_values collection
    if (length(input_files) > 0) {
      if ("--input" %in% names(array_values)) {
        array_values[["--input"]] <- c(array_values[["--input"]], input_files)
      } else {
        array_values[["--input"]] <- input_files
      }
    }
    
    if (length(sample_names) > 0) {
      if ("--name" %in% names(array_values)) {
        array_values[["--name"]] <- c(array_values[["--name"]], sample_names)
      } else {
        array_values[["--name"]] <- sample_names
      }
    }
    
    message("Reclassifying positional arguments:")
    message("  Input files:", paste(input_files, collapse=", "))
    message("  Sample names:", paste(sample_names, collapse=", "))
  }
  
  # Build new args list
  new_args <- c()
  
  # Add regular args
  for (arg_name in names(other_args)) {
    new_args <- c(new_args, arg_name, other_args[[arg_name]])
  }
  
  # Add required args if found
  for (req_arg in required_args) {
    if (req_arg %in% names(required_values)) {
      new_args <- c(new_args, req_arg, required_values[[req_arg]])
    }
  }
  
  # Add array args
  for (arg_name in names(array_values)) {
    for (value in array_values[[arg_name]]) {
      new_args <- c(new_args, arg_name, value)
    }
  }
  
  # Check if we have all the key arguments
  if (length(new_args) > 0) {
    message("Reformatted argument list:")
    i <- 1
    while (i <= length(new_args)) {
      if (i+1 <= length(new_args)) {
        message(paste("  ", new_args[i], new_args[i+1]))
        i <- i + 2
      } else {
        message(paste("  ", new_args[i], "TRUE"))
        i <- i + 1
      }
    }
    return(new_args)
  }
  
  # Return original arguments if no preprocessing needed
  return(all_args)
}

# Source utility functions from common directory
# First try Docker standard path, then fall back to relative path
if (file.exists("/usr/local/bin/functions/common/utilities.R")) {
  source("/usr/local/bin/functions/common/utilities.R")
} else if (file.exists("functions/common/utilities.R")) {
  source("functions/common/utilities.R")
} else {
  stop("Could not find utilities.R file")
}

# Source all required function files
report_memory_usage("Before loading function files")

# Source common functions
source_with_fallback("functions/common/visualization.R", "/usr/local/bin/functions/common/visualization.R")
source_with_fallback("functions/common/export_functions.R", "/usr/local/bin/functions/common/export_functions.R")

# Source DESeq2 LRT Step 1 specific functions
source_with_fallback("functions/deseq2_lrt_step_1/cli_args.R", "/usr/local/bin/functions/deseq2_lrt_step_1/cli_args.R")
source_with_fallback("functions/deseq2_lrt_step_1/data_processing.R", "/usr/local/bin/functions/deseq2_lrt_step_1/data_processing.R")
source_with_fallback("functions/deseq2_lrt_step_1/deseq2_analysis.R", "/usr/local/bin/functions/deseq2_lrt_step_1/deseq2_analysis.R")
source_with_fallback("functions/deseq2_lrt_step_1/contrast_generation.R", "/usr/local/bin/functions/deseq2_lrt_step_1/contrast_generation.R")
source_with_fallback("functions/deseq2_lrt_step_1/workflow.R", "/usr/local/bin/functions/deseq2_lrt_step_1/workflow.R")

report_memory_usage("After loading function files")

# Configure plot theme
configure_plot_theme()

# Main execution with memory management and error handling
main_with_memory_management <- function() {
  tryCatch({
    # Pre-process command line arguments first - improves compatibility with CWL
    processed_args <- preprocess_args()

    # Override command line arguments with processed ones
    # This is a hack to make the argparse parser work with our processed arguments
    if (!identical(processed_args, commandArgs(trailingOnly = TRUE))) {
      message("Command line arguments have been preprocessed for compatibility")
      
      # Create a safer version of the command line arguments
      safe_args <- c()
      for (arg in processed_args) {
        # Escape any quotes and special characters
        escaped_arg <- gsub('"', '\\\\"', arg)
        # Add quotes around arguments with spaces
        if (grepl(" ", escaped_arg)) {
          escaped_arg <- paste0('"', escaped_arg, '"')
        }
        safe_args <- c(safe_args, escaped_arg)
      }
      
      # Temporarily save the processed arguments to a file and source it
      temp_args_file <- tempfile(pattern = "args_", fileext = ".R")
      on.exit(unlink(temp_args_file), add = TRUE)
      
      cat("commandArgs_original <- commandArgs\n", file = temp_args_file)
      cat("commandArgs <- function(trailingOnly = FALSE) {\n", file = temp_args_file, append = TRUE)
      cat("  if (trailingOnly) {\n", file = temp_args_file, append = TRUE)
      
      # Convert the arguments to a proper R character vector
      arg_vector <- paste("c(", paste(shQuote(processed_args), collapse = ", "), ")")
      cat(paste0("    return(", arg_vector, ")\n"), file = temp_args_file, append = TRUE)
      
      cat("  } else {\n", file = temp_args_file, append = TRUE)
      cat("    return(commandArgs_original(FALSE))\n", file = temp_args_file, append = TRUE)
      cat("  }\n", file = temp_args_file, append = TRUE)
      cat("}\n", file = temp_args_file, append = TRUE)
      
      source(temp_args_file)
    }

    ### Parse CLI arguments
    message("Parsing command line arguments...")
    args <- tryCatch({
      params::get_cli_args()
    }, error = function(e) {
      message("Error during argument parsing: ", e$message)
      
      # Emergency fallback for required arguments
      message("EMERGENCY: Attempting to extract required arguments")
      
      # Get original command line args
      orig_args <- commandArgs(trailingOnly = TRUE)
      
      # Initialize empty result with NULL values
      emergency_args <- list(
        meta = NULL,
        design = NULL,
        reduced = NULL,
        input = character(0),
        name = character(0)
      )
      
      # Extract required arguments directly if parser failed
      # This is a more aggressive approach to handle malformed arguments
      required_flags <- c("--meta", "--design", "--reduced", "--input", "--name")
      
      # Handle regular args that aren't array-based
      for (flag in c("--meta", "--design", "--reduced")) {
        flag_idx <- which(orig_args == flag)
        if (length(flag_idx) > 0 && flag_idx[1] < length(orig_args)) {
          value <- orig_args[flag_idx[1] + 1]
          if (!grepl("^--", value)) {
            param_name <- gsub("^--", "", flag)
            emergency_args[[param_name]] <- value
            message("Found ", flag, " = ", value)
          }
        }
      }
      
      # Handle array args
      for (flag in c("--input", "--name")) {
        flag_indices <- which(orig_args == flag)
        values <- character(0)
        
        # Get all values for the flag
        for (idx in flag_indices) {
          if (idx < length(orig_args) && !grepl("^--", orig_args[idx + 1])) {
            values <- c(values, orig_args[idx + 1])
          }
        }
        
        # Store values if found
        if (length(values) > 0) {
          param_name <- gsub("^--", "", flag)
          emergency_args[[param_name]] <- values
          message("Found ", flag, " = ", paste(values, collapse=", "))
        }
      }
      
      # Check for positional arguments that could be inputs or sample names
      positional_args <- orig_args[!grepl("^--", orig_args)]
      positional_args <- setdiff(positional_args, unlist(emergency_args))
      
      if (length(positional_args) > 0) {
        message("Found ", length(positional_args), " potential positional arguments")
        
        # Try to identify input files vs sample names
        file_pattern <- "\\.(tsv|csv)$"
        input_files <- positional_args[grepl(file_pattern, positional_args)]
        sample_names <- positional_args[!grepl(file_pattern, positional_args)]
        
        if (length(input_files) > 0) {
          emergency_args$input <- c(emergency_args$input, input_files)
          message("Adding ", length(input_files), " input files from positional args")
        }
        
        if (length(sample_names) > 0) {
          emergency_args$name <- c(emergency_args$name, sample_names)
          message("Adding ", length(sample_names), " sample names from positional args")
        }
      }
      
      # Check if we found all required arguments
      missing_args <- names(emergency_args)[sapply(emergency_args, function(x) is.null(x) || (is.character(x) && length(x) == 0))]
      if (length(missing_args) > 0) {
        message("ERROR: Still missing required arguments: ", paste(missing_args, collapse=", "))
        stop("Missing required arguments: ", paste(missing_args, collapse=", "))
      }
      
      # Set defaults for other parameters
      default_params <- list(
        output = "results",
        rpkm_cutoff = 1,
        k_hopach = 5,
        kmax_hopach = 6,
        cosine_normalized = TRUE,
        distance_method = "cosine",
        clustering_method = "complete",
        filter = TRUE,
        batch_correct = FALSE,
        center_rowmeans = TRUE,
        fdr = 0.1
      )
      
      # Merge with existing args
      for (param in names(default_params)) {
        if (is.null(emergency_args[[param]])) {
          emergency_args[[param]] <- default_params[[param]]
        }
      }
      
      # Print final parameters
      message("Emergency argument parsing completed!")
      message("Parameters:")
      for (param in names(emergency_args)) {
        param_value <- emergency_args[[param]]
        if (is.character(param_value) && length(param_value) > 1) {
          message(paste0("  ", param, ": [", paste(param_value, collapse=", "), "]"))
        } else {
          message(paste0("  ", param, ": ", param_value))
        }
      }
      
      return(emergency_args)
    })
    
    ### Create output dir if not exist
    if (!dir.exists(args$output)) {
      dir.create(args$output, recursive = TRUE)
    }
    
    ### Setup logging
    logfile_path <- file.path(args$output, "log.txt")
    logger::log_appender(logger::appender_tee(logfile_path))
    set_log_level_verbose()
    
    log_info("Starting DESeq2 analysis with LRT")
    log_debug("Arguments:", print_all_args(args))
    
    ### Validate CLI arguments
    params::assert_args(args)
    log_debug("CLI arguments validated")
    
    ### Run main process
    run_main_process(args)
  }, error = function(e) {
    # Also log to stderr - important for debugging CWL
    write(paste("ERROR:", e$message), stderr())
    
    # Try to log to file if possible
    tryCatch({
      log_error(e$message)
      log_error(e$call)
      log_error(traceback())
    }, error = function(e) {
      write(paste("Failed to log error to file:", e$message), stderr())
    })
    
    # Make sure there's a non-zero exit 
    quit(status = 1)
  })
}

# Execute main function
main_with_memory_management()
