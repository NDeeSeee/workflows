#!/usr/bin/env Rscript
#
# DESeq2 Differential Expression Analysis
#
# This script performs differential expression analysis using DESeq2.
# It's been refactored for better maintainability with functions organized
# into separate files.
#
# Version: 0.1.0

# First load the utilities module which contains the source_with_fallback function
if (file.exists("functions/common/utilities.R")) {
  source("functions/common/utilities.R")
} else if (file.exists("/usr/local/bin/functions/common/utilities.R")) {
  source("/usr/local/bin/functions/common/utilities.R")
} else {
  # If we can't find utilities, define source_with_fallback function minimally to continue
  source_with_fallback <- function(file_path, fallback_path) {
    if (file.exists(file_path)) {
      source(file_path)
      return(TRUE)
    } else if (!is.null(fallback_path) && file.exists(fallback_path)) {
      source(fallback_path)
      return(TRUE)
    }
    return(FALSE)
  }
}

# Source the workflow file
if (!source_with_fallback("functions/deseq/workflow.R", "/usr/local/bin/functions/deseq/workflow.R")) {
  stop("Could not find workflow.R file")
}

# Initialize the environment (loads libraries, sets options, loads source files)
initialize_environment()

# Execute main function
main_with_memory_management()
