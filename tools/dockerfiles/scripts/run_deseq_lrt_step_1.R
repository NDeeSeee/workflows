#!/usr/bin/env Rscript
#
# DESeq2 LRT Analysis - Step 1
#
# This script performs differential expression analysis using DESeq2 with 
# Likelihood Ratio Test (LRT). It's been refactored for better maintainability
# with functions organized into separate files.
#
# Version: 0.1.5

# Define a basic error handler
options(error = function() {
  message("\nExecution halted due to error: ", geterrmessage())
  quit(status = 1)
})

# Set up a basic log message function in case the real one isn't available
if (!exists("log_message")) {
  log_message <- function(message, level="INFO") {
    cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), level, ":", message, "\n")
  }
}

# First load the utilities module which contains the source_with_fallback function
tryCatch({
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
  result <- source_with_fallback("functions/deseq2_lrt_step_1/workflow.R", "/usr/local/bin/functions/deseq2_lrt_step_1/workflow.R")
  if (!result) {
    stop("Could not find workflow.R file")
  }

  # Initialize the environment (loads libraries, sets options, loads source files)
  initialize_environment()

  # Execute main function
  main_with_memory_management()
}, error = function(e) {
  message("Error during initialization: ", e$message)
  quit(status = 1)
})
