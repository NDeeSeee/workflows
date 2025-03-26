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
  
  # Check if there are any arguments without -- prefix
  has_positional <- any(!grepl("^--", all_args))
  
  if (has_positional) {
    # Extract positional arguments
    positional_args <- all_args[!grepl("^--", all_args) & !grepl("^-", all_args)]
    
    # Separate positional arguments into names and inputs based on file path pattern
    file_pattern <- "\\.(tsv|csv)$"
    input_files <- positional_args[grepl(file_pattern, positional_args)]
    sample_names <- positional_args[!grepl(file_pattern, positional_args)]
    
    # If we have both names and inputs, filter out these positional arguments and add them with proper flags
    if (length(input_files) > 0 && length(sample_names) > 0) {
      # Remove all positional arguments
      flag_args <- all_args[grepl("^--", all_args) | grepl("^-", all_args)]
      flag_values <- all_args[which(grepl("^--", all_args) | grepl("^-", all_args)) + 1]
      flag_values <- flag_values[!grepl("^--", flag_values) & !grepl("^-", flag_values)]
      
      # Combine flags with their values
      flag_pairs <- c()
      flag_index <- 1
      for (i in 1:length(flag_args)) {
        flag_pairs <- c(flag_pairs, flag_args[i])
        if (flag_index <= length(flag_values) && 
            !grepl("^--", flag_values[flag_index]) && 
            !grepl("^-", flag_values[flag_index])) {
          flag_pairs <- c(flag_pairs, flag_values[flag_index])
          flag_index <- flag_index + 1
        }
      }
      
      # Create new argument vector with proper flags
      new_args <- flag_pairs
      
      # Add input files with proper flags
      for (file in input_files) {
        new_args <- c(new_args, "--input", file)
      }
      
      # Add sample names with proper flags
      for (name in sample_names) {
        new_args <- c(new_args, "--name", name)
      }
      
      # Replace command line arguments
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
    
    # Get and parse command line arguments
    args <- get_args()
    
    # Run the main analysis pipeline
    run_deseq_analysis(args)
    
    log_message("DESeq2 LRT Step 1 completed successfully", "SUCCESS")
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
