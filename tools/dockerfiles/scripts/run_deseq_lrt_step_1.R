#!/usr/bin/env Rscript
#
# DESeq2 LRT Analysis - Step 1
#
# This script performs differential expression analysis using DESeq2 with 
# Likelihood Ratio Test (LRT). It's been refactored for better maintainability
# with functions organized into separate files.
#
# Version: 0.1.5

# Load the workflow
source("functions/deseq2_lrt_step_1/workflow.R")

# Initialize the environment (loads libraries, sets options, loads source files)
initialize_environment()

# Execute main function
main_with_memory_management()
