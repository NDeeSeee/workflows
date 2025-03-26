#!/usr/bin/env Rscript
#
# DESeq2 LRT Analysis - Step 2
#
# This script performs contrast-specific analyses using the DESeq2 object created in Step 1.
# It processes selected contrasts, generates visualizations, and exports results in various formats.
#
# Version: 0.1.0

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
  
  # For clustering
  library(hopach)
  
  # For visualization
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(ggplot2)
  library(ggrepel)
  library(plotly)
  
  # For GCT export
  library(cmapR)
  
  # For memory profiling
  library(pryr)
})

# Use dplyr functions with proper namespace to avoid conflicts
mutate <- dplyr::mutate
filter <- dplyr::filter
group_by <- dplyr::group_by
slice <- dplyr::slice
rename <- dplyr::rename
select <- dplyr::select
arrange <- dplyr::arrange
distinct <- dplyr::distinct
`%>%` <- magrittr::`%>%`
`%in%` <- base::`%in%`
`%/%` <- base::`%/%`
`%%` <- base::`%%`

# Source the utilities file first (needed for source_with_fallback function)
if (file.exists("/usr/local/bin/functions/common/utilities.R")) {
  source("/usr/local/bin/functions/common/utilities.R")
} else if (file.exists("functions/common/utilities.R")) {
  source("functions/common/utilities.R")
} else {
  stop("Could not find utilities.R file")
}

# Now we can use report_memory_usage and source_with_fallback from utilities.R
report_memory_usage("Before loading common functions")

# Source common additional function modules
source_with_fallback("functions/common/constants.R", "/usr/local/bin/functions/common/constants.R")
source_with_fallback("functions/common/error_handling.R", "/usr/local/bin/functions/common/error_handling.R")
source_with_fallback("functions/common/logging.R", "/usr/local/bin/functions/common/logging.R")
source_with_fallback("functions/common/visualization.R", "/usr/local/bin/functions/common/visualization.R")
source_with_fallback("functions/common/clustering.R", "/usr/local/bin/functions/common/clustering.R")
source_with_fallback("functions/common/export_functions.R", "/usr/local/bin/functions/common/export_functions.R")
report_memory_usage("After loading common functions")

# Source DESeq2 LRT Step 2 specific functions
source_with_fallback("functions/deseq2_lrt_step_2/cli_args.R", "/usr/local/bin/functions/deseq2_lrt_step_2/cli_args.R")
source_with_fallback("functions/deseq2_lrt_step_2/data_processing.R", "/usr/local/bin/functions/deseq2_lrt_step_2/data_processing.R")
source_with_fallback("functions/deseq2_lrt_step_2/contrast_analysis.R", "/usr/local/bin/functions/deseq2_lrt_step_2/contrast_analysis.R")
source_with_fallback("functions/deseq2_lrt_step_2/workflow.R", "/usr/local/bin/functions/deseq2_lrt_step_2/workflow.R")

# Configure plot theme
configure_plot_theme()

# Wrapper function with memory management
main_with_memory_management <- function() {
  # Start timing
  start_time <- Sys.time()
  log_message(paste("DESeq2 LRT Step 2 started at", format(start_time, "%Y-%m-%d %H:%M:%S")))
  
  # Get command line arguments
  args <- with_error_handling({
    get_args()
  })
  
  if (is.null(args)) {
    log_error("Failed to parse command line arguments")
    stop("Failed to parse command line arguments")
  }
  
  # Configure parallel processing
  if (args$threads > 1) {
    log_message(paste("Setting up parallel execution with", args$threads, "threads"))
    register(MulticoreParam(args$threads))
  } else {
    log_message("Running in single-threaded mode")
  }
  
  # Run the main workflow
  results <- with_error_handling({
    run_workflow(args)
  })
  
  if (is.null(results)) {
    log_error("Workflow execution failed")
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  # Report end time and duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  log_message(paste("DESeq2 LRT Step 2 completed at", format(end_time, "%Y-%m-%d %H:%M:%S")))
  log_message(paste("Total execution time:", round(as.numeric(duration), 2), "minutes"))
  
  # Clean up large objects to free memory
  rm(results)
  invisible(gc())
  
  # Final memory report
  report_memory_usage("Final")
  
  log_message("DESeq2 LRT Step 2 analysis completed successfully")
}

# Execute main function with memory management
main_with_memory_management()
