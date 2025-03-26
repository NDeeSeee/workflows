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
