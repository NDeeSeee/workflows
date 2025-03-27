#!/usr/bin/env Rscript
#
# DESeq2 LRT Analysis - Step 2
#
# This script performs contrast-specific analyses using the DESeq2 object created in Step 1.
# It processes selected contrasts, generates visualizations, and exports results in various formats.
#
# Version: 0.1.0

source("functions/deseq2_lrt_step_2/workflow.R")

# Initialize the environment (loads libraries, sets options, loads source files)
initialize_environment()

# Execute main function
main_with_memory_management()
