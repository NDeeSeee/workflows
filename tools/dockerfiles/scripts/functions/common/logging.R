#!/usr/bin/env Rscript

# --- Logging functions ---

# Debug mode control
DEBUG_MODE <- FALSE

# Function to enable debug mode
enable_debug_mode <- function(args) {
  assign("DEBUG_MODE", TRUE, envir = .GlobalEnv)
  log_message("Debug mode enabled", "DEBUG")
  
  if (!is.null(args$verbose) && args$verbose) {
    log_message("Verbose logging enabled", "DEBUG")
    # Configure more detailed error messages for debugging
    options(error = function() {
      traceback(2)
      message("\nError: ", geterrmessage())
      q(status = 1)
    })
  }
  
  # Clean up any existing log file at the start
  if (args$clean_logs) {
    file.remove("deseq_analysis.log")
    log_message("Started new log file", "INFO")
  }
}

# Formatted logging function with log levels
log_message <- function(message, level = "INFO") {
  # Define colors for different log levels
  colors <- list(
    "DEBUG" = "\033[36m",  # Cyan
    "INFO" = "\033[0m",    # Default
    "STEP" = "\033[34m",   # Blue
    "WARNING" = "\033[33m", # Yellow
    "ERROR" = "\033[31m",  # Red
    "SUCCESS" = "\033[32m" # Green
  )
  
  # Reset color
  reset <- "\033[0m"
  
  # Format timestamp
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Create formatted message
  log_string <- paste0(
    "[", timestamp, "] ", 
    "[", level, "] ",
    message
  )
  
  # Print to console with color
  formatted_message <- paste0(
    colors[[level]], 
    log_string,
    reset
  )
  cat(formatted_message, "\n")
  
  # Also append to log file
  write(log_string, file = "deseq_analysis.log", append = TRUE)
}

# Debug-specific logging that only appears when debug mode is on
debug_log <- function(message, data = NULL) {
  if (exists("DEBUG_MODE") && DEBUG_MODE) {
    log_message(message, "DEBUG")
    
    # If data object is provided, print its summary
    if (!is.null(data)) {
      if (is.data.frame(data)) {
        log_message(glue::glue("Data dimensions: {nrow(data)} rows x {ncol(data)} columns"), "DEBUG")
        if (nrow(data) > 0) {
          class_info <- sapply(data, class)
          log_message("Column types: ", "DEBUG")
          print(class_info)
          
          # Print a few rows for inspection
          log_message("Preview of data:", "DEBUG")
          print(head(data, 5))
        }
      } else if (is.matrix(data)) {
        log_message(glue::glue("Matrix dimensions: {nrow(data)} rows x {ncol(data)} columns"), "DEBUG")
        log_message("Matrix preview:", "DEBUG")
        print(head(data, 5))
      } else if (is.list(data)) {
        log_message(glue::glue("List with {length(data)} elements"), "DEBUG")
        log_message(glue::glue("List names: {paste(names(data), collapse=', ')}"), "DEBUG")
      } else {
        log_message("Object summary:", "DEBUG")
        print(summary(data))
      }
    }
  }
} 