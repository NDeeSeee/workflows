#!/usr/bin/env Rscript
#
# DESeq2 LRT Analysis - Step 1
#
# This script performs differential expression analysis using DESeq2 with 
# Likelihood Ratio Test (LRT). It's been refactored for better maintainability
# with functions organized into separate files.
#
# Version: 0.1.1

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
  # For argument parsing
  library(argparse)
  
  # For parallel processing
  library(BiocParallel)
  
  # For DESeq2 analysis
  library(DESeq2)
  
  # For data manipulation
  library(tidyverse)
  library(data.table)
  
  # For batch correction
  library(limma)
  
  # For visualization
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(ggplot2)
  library(ggrepel)
  library(plotly)
  
  # For GCT export
  library(cmapR)
  
  # For utilities
  library(rlang)
  library(stringr)
  library(glue)
})

# Use dplyr functions with proper namespace to avoid conflicts
mutate <- dplyr::mutate
filter <- dplyr::filter
group_by <- dplyr::group_by
slice <- dplyr::slice
rename <- dplyr::rename
select <- dplyr::select
`%>%` <- magrittr::`%>%`
`%in%` <- base::`%in%`

# Memory usage tracking function
report_memory_usage <- function(label = "") {
  gc(verbose = FALSE)
  mem_used <- pryr::mem_used()
  message(paste0("[Memory] ", label, ": ", round(mem_used / 1024^2, 1), " MB"))
}

# Helper function to source files with fallback paths
source_with_fallback <- function(filepath, absolute_path = NULL) {
  # Try absolute path first if provided
  if (!is.null(absolute_path) && file.exists(absolute_path)) {
    message(paste("Sourcing from absolute path:", absolute_path))
    return(source(absolute_path))
  }
  
  # Try relative path
  if (file.exists(filepath)) {
    message(paste("Sourcing from relative path:", filepath))
    return(source(filepath))
  }
  
  # If we get here, try a standard Docker path
  docker_path <- file.path("/usr/local/bin", filepath)
  if (file.exists(docker_path)) {
    message(paste("Sourcing from Docker path:", docker_path))
    return(source(docker_path))
  }
  
  # If all fails, error
  stop(paste("Could not find file to source:", filepath))
}

# Source utility functions
report_memory_usage("Before loading common utilities")
source_with_fallback("functions/common/utilities.R", "/usr/local/bin/functions/common/utilities.R")
report_memory_usage("After loading common utilities")

# Source common visualization and export functions
source_with_fallback("functions/common/visualization.R", "/usr/local/bin/functions/common/visualization.R")
source_with_fallback("functions/common/export_functions.R", "/usr/local/bin/functions/common/export_functions.R")
report_memory_usage("After loading common visualization and export functions")

# Source DESeq2 LRT Step 1 specific functions
source_with_fallback("functions/deseq2_lrt_step_1/cli_args.R", "/usr/local/bin/functions/deseq2_lrt_step_1/cli_args.R")
source_with_fallback("functions/deseq2_lrt_step_1/data_processing.R", "/usr/local/bin/functions/deseq2_lrt_step_1/data_processing.R")
source_with_fallback("functions/deseq2_lrt_step_1/deseq2_analysis.R", "/usr/local/bin/functions/deseq2_lrt_step_1/deseq2_analysis.R")
source_with_fallback("functions/deseq2_lrt_step_1/contrast_generation.R", "/usr/local/bin/functions/deseq2_lrt_step_1/contrast_generation.R")
source_with_fallback("functions/deseq2_lrt_step_1/workflow.R", "/usr/local/bin/functions/deseq2_lrt_step_1/workflow.R")

# Configure plot theme
configure_plot_theme()

# Wrapper function with memory management
main_with_memory_management <- function() {
  # Start timing
  start_time <- Sys.time()
  message(glue::glue("DESeq2 LRT Step 1 started at {format(start_time, '%Y-%m-%d %H:%M:%S')}"))
  
  # Get command line arguments
  args <- get_args()
  
  # Configure parallel processing
  if (args$threads > 1) {
    message(paste("Setting up parallel execution with", args$threads, "threads"))
    register(MulticoreParam(args$threads))
  } else {
    message("Running in single-threaded mode")
  }
  
  # Run the main workflow with validated args
  main(args)
  
  # Report end time and duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "secs")
  message(glue::glue("DESeq2 LRT Step 1 completed at {format(end_time, '%Y-%m-%d %H:%M:%S')}"))
  message(glue::glue("Total execution time: {round(as.numeric(duration), 2)} seconds"))
  
  # Final memory report
  report_memory_usage("Final")
}

# Execute main function with memory management
main_with_memory_management()
