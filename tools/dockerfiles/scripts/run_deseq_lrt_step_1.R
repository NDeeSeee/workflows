#!/usr/bin/env Rscript
#
# DESeq2 LRT Analysis - Step 1
#
# This script performs differential expression analysis using DESeq2 with 
# Likelihood Ratio Test (LRT). It's been refactored for better maintainability
# with functions organized into separate files.
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
  
  # For batch correction
  library(sva)
  
  # For visualization
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
  library(Glimma)
  
  # For data export
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

# Source function files from the functions directory
# Common functions
source("/usr/local/bin/functions/common/logging.R")
source("/usr/local/bin/functions/common/error_handling.R")
source("/usr/local/bin/functions/common/utilities.R")

# DESeq2 LRT Step 1 specific functions
source("/usr/local/bin/functions/deseq2_lrt_step_1/cli_args.R")
source("/usr/local/bin/functions/deseq2_lrt_step_1/data_processing.R")
source("/usr/local/bin/functions/deseq2_lrt_step_1/metadata_validation.R")
source("/usr/local/bin/functions/deseq2_lrt_step_1/deseq2_analysis.R")
source("/usr/local/bin/functions/deseq2_lrt_step_1/contrast_generation.R")
source("/usr/local/bin/functions/deseq2_lrt_step_1/export_functions.R")
source("/usr/local/bin/functions/deseq2_lrt_step_1/visualization.R")
source("/usr/local/bin/functions/deseq2_lrt_step_1/workflow.R")

# Configure plot theme
configure_plot_theme()

# Execute main function
main()
