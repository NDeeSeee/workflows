#!/usr/bin/env Rscript
#
# DESeq2 Differential Expression Analysis
#
# This script performs differential expression analysis using DESeq2.
# It's been refactored for better maintainability with functions organized
# into separate files.
#
# Version: 0.1.0


# Load the workflow
source("functions/deseq/workflow.R")

# Initialize the environment (loads libraries, sets options, loads source files)
initialize_environment()

# Execute main function
main_with_memory_management()
