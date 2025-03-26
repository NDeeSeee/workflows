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
  
  # Extract key required arguments first (meta, design, reduced)
  required_args <- c("--meta", "--design", "--reduced")
  required_values <- list()
  for (req_arg in required_args) {
    arg_index <- which(all_args == req_arg)
    if (length(arg_index) > 0 && arg_index[1] < length(all_args)) {
      required_values[[req_arg]] <- all_args[arg_index[1] + 1]
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
    
    # If we have both names and inputs as positional args, create a proper argument list
    if (length(input_files) > 0 || length(sample_names) > 0) {
      message("Reclassifying positional arguments:")
      message("  Input files:", paste(input_files, collapse=", "))
      message("  Sample names:", paste(sample_names, collapse=", "))
      
      # Start with non-positional arguments
      new_args <- c()
      
      # Add all flag args and their values
      i <- 1
      while (i <= length(all_args)) {
        if (grepl("^--", all_args[i]) || grepl("^-[a-zA-Z]", all_args[i])) {
          new_args <- c(new_args, all_args[i])
          
          # If next arg exists and is not a flag, it's a value - add it too
          if (i < length(all_args) && 
              !grepl("^--", all_args[i+1]) && 
              !grepl("^-[a-zA-Z]", all_args[i+1])) {
            new_args <- c(new_args, all_args[i+1])
            i <- i + 2
          } else {
            # Just a flag without value
            i <- i + 1
          }
        } else {
          # Skip positional args here - we'll add them properly later
          i <- i + 1
        }
      }
      
      # Now add properly flagged input files
      for (file in input_files) {
        new_args <- c(new_args, "--input", file)
      }
      
      # Now add properly flagged sample names
      for (name in sample_names) {
        new_args <- c(new_args, "--name", name)
      }
      
      # Check if we have the required arguments
      for (req_arg in required_args) {
        if (!req_arg %in% new_args && !is.null(required_values[[req_arg]])) {
          # Add the required arg and its value if we found it earlier
          new_args <- c(new_args, req_arg, required_values[[req_arg]])
        }
      }
      
      message("Reformatted argument list:")
      for (i in seq(1, length(new_args), by=2)) {
        if (i+1 <= length(new_args)) {
          message(paste("  ", new_args[i], new_args[i+1]))
        } else {
          message(paste("  ", new_args[i], "TRUE"))
        }
      }
      
      return(new_args)
    }
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
  log_message("DESeq2 LRT Step 1 started", "START")
  
  tryCatch({
    # Pre-process command line arguments to handle unexpected formats
    processed_args <- preprocess_args()
    
    # Override command line arguments with processed ones
    # This is a hack to make the argparse parser work with our processed arguments
    if (!identical(processed_args, commandArgs(trailingOnly = TRUE))) {
      message("Command line arguments have been preprocessed for compatibility")
      
      # Temporarily save the processed arguments to a file and source it
      temp_args_file <- tempfile(pattern = "args_", fileext = ".R")
      on.exit(unlink(temp_args_file), add = TRUE)
      
      cat("commandArgs_original <- commandArgs\n", file = temp_args_file)
      cat("commandArgs <- function(trailingOnly = FALSE) {\n", file = temp_args_file, append = TRUE)
      cat("  if (trailingOnly) {\n", file = temp_args_file, append = TRUE)
      cat("    return(c(", paste0("\"", processed_args, "\"", collapse = ", "), "))\n", file = temp_args_file, append = TRUE)
      cat("  } else {\n", file = temp_args_file, append = TRUE)
      cat("    return(commandArgs_original(FALSE))\n", file = temp_args_file, append = TRUE)
      cat("  }\n", file = temp_args_file, append = TRUE)
      cat("}\n", file = temp_args_file, append = TRUE)
      
      source(temp_args_file)
    }
    
    # Try to parse command line arguments
    tryCatch({
      # Get and parse command line arguments
      args <- get_args()
      
      # Additional emergency checks for required arguments
      required_args <- c("meta", "design", "reduced", "input", "name")
      missing_args <- required_args[!required_args %in% names(args)]
      
      if (length(missing_args) > 0) {
        message(paste("ERROR: Still missing required arguments after processing:", paste(missing_args, collapse=", ")))
        
        # Emergency override if we're in a really bad state
        if ("meta" %in% missing_args || "design" %in% missing_args || "reduced" %in% missing_args) {
          # Try to extract these from the original arguments as a last resort
          all_args <- commandArgs(trailingOnly = TRUE)
          
          # Look for --meta
          meta_idx <- which(all_args == "--meta")
          if (length(meta_idx) > 0 && meta_idx[1] < length(all_args)) {
            args$meta <- all_args[meta_idx[1] + 1]
            message(paste("Emergency override: Setting meta =", args$meta))
          }
          
          # Look for --design
          design_idx <- which(all_args == "--design")
          if (length(design_idx) > 0 && design_idx[1] < length(all_args)) {
            args$design <- all_args[design_idx[1] + 1]
            message(paste("Emergency override: Setting design =", args$design))
          }
          
          # Look for --reduced
          reduced_idx <- which(all_args == "--reduced")
          if (length(reduced_idx) > 0 && reduced_idx[1] < length(all_args)) {
            args$reduced <- all_args[reduced_idx[1] + 1]
            message(paste("Emergency override: Setting reduced =", args$reduced))
          }
        }
        
        # Input files and sample names are trickier, look for positional args
        if ("input" %in% missing_args || "name" %in% missing_args) {
          all_args <- commandArgs(trailingOnly = TRUE)
          
          # Look for file patterns
          file_pattern <- "\\.(tsv|csv)$"
          potential_files <- all_args[grepl(file_pattern, all_args)]
          
          # Look for potential sample names
          potential_names <- all_args[!grepl("^--", all_args) & !grepl(file_pattern, all_args)]
          
          # Set input and name if we found them
          if (length(potential_files) > 0) {
            args$input <- potential_files
            message(paste("Emergency override: Setting input to", length(potential_files), "files"))
          }
          
          if (length(potential_names) > 0) {
            args$name <- potential_names
            message(paste("Emergency override: Setting name to", length(potential_names), "names"))
          }
        }
        
        # One final check
        still_missing <- required_args[!required_args %in% names(args)]
        if (length(still_missing) > 0) {
          stop(paste("Cannot continue - still missing required arguments:", paste(still_missing, collapse=", ")))
        }
      }
      
      # Run the main analysis pipeline
      run_deseq_analysis(args)
      
      log_message("DESeq2 LRT Step 1 completed successfully", "SUCCESS")
    }, error = function(e) {
      # Special handling for argparse errors
      if (grepl("parse error", e$message) || grepl("unrecognized arguments", e$message)) {
        message("ERROR: Failed to parse arguments automatically. Trying manual fallback...")
        
        # Manual argument extraction from command line
        all_args <- commandArgs(trailingOnly = TRUE)
        
        # Create an empty args list
        args <- list()
        
        # Extract required arguments directly
        meta_idx <- which(all_args == "--meta")
        if (length(meta_idx) > 0 && meta_idx[1] < length(all_args)) {
          args$meta <- all_args[meta_idx[1] + 1]
        }
        
        design_idx <- which(all_args == "--design")
        if (length(design_idx) > 0 && design_idx[1] < length(all_args)) {
          args$design <- all_args[design_idx[1] + 1]
        }
        
        reduced_idx <- which(all_args == "--reduced")
        if (length(reduced_idx) > 0 && reduced_idx[1] < length(all_args)) {
          args$reduced <- all_args[reduced_idx[1] + 1]
        }
        
        # Handle input files and sample names
        # Look for file patterns
        file_pattern <- "\\.(tsv|csv)$"
        args$input <- all_args[grepl(file_pattern, all_args)]
        
        # Find all non-flag, non-file arguments that might be sample names
        potential_names <- all_args[!grepl("^--", all_args) & !grepl("^-", all_args) & !grepl(file_pattern, all_args)]
        
        # Filter out argument values
        flag_indices <- which(grepl("^--", all_args))
        arg_values <- c()
        for (idx in flag_indices) {
          if (idx < length(all_args) && !grepl("^--", all_args[idx+1])) {
            arg_values <- c(arg_values, all_args[idx+1])
          }
        }
        
        args$name <- setdiff(potential_names, arg_values)
        
        # Extract other common parameters
        batchcorrection_idx <- which(all_args == "--batchcorrection")
        if (length(batchcorrection_idx) > 0 && batchcorrection_idx[1] < length(all_args)) {
          args$batchcorrection <- all_args[batchcorrection_idx[1] + 1]
        } else {
          args$batchcorrection <- "none"
        }
        
        output_idx <- which(all_args == "--output")
        if (length(output_idx) > 0 && output_idx[1] < length(all_args)) {
          args$output_prefix <- all_args[output_idx[1] + 1]
        } else {
          args$output_prefix <- "./deseq_lrt_step_1"
        }
        
        # Check if we have all required arguments
        required_args <- c("meta", "design", "reduced", "input", "name")
        missing_args <- required_args[!required_args %in% names(args)]
        
        if (length(missing_args) > 0) {
          stop(paste("Cannot continue - missing required arguments:", paste(missing_args, collapse=", ")))
        }
        
        # Set defaults for other parameters
        args$fdr <- 0.1
        args$lfcthreshold <- 0.59
        args$cluster_method <- "none"
        args$threads <- 1
        
        # Run the analysis with manually extracted args
        message("Running with manually extracted arguments:")
        for (key in names(args)) {
          if (length(args[[key]]) > 1) {
            message(paste("  ", key, ": [array of", length(args[[key]]), "items]"))
          } else {
            message(paste("  ", key, ":", args[[key]]))
          }
        }
        
        run_deseq_analysis(args)
        log_message("DESeq2 LRT Step 1 completed successfully (with manual argument handling)", "SUCCESS")
      } else {
        # Rethrow the error for other types of errors
        stop(e)
      }
    })
  }, error = function(e) {
    log_message(paste("ERROR:", conditionMessage(e)), "ERROR")
    log_message(paste("See traceback for details:", deparse(e$call)), "ERROR")
    print(rlang::last_trace())
    # Exit with error code
    quit(save = "no", status = 1, runLast = FALSE)
  }, finally = {
    # Clean up
    gc()
  })
}

# Execute main function
main_with_memory_management()
