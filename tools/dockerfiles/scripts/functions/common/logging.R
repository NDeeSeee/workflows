#!/usr/bin/env Rscript
#
# Common logging functions
#

#' Log a message to the console
#'
#' @param message The message to log
#' @param timestamp Whether to include a timestamp
#' @export
log_message <- function(message, timestamp = TRUE) {
  if (timestamp) {
    timestamp_str <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    cat(timestamp_str, "INFO:", message, "\n")
  } else {
    cat("INFO:", message, "\n")
  }
}

#' Log a warning message to the console
#'
#' @param message The warning message to log
#' @param timestamp Whether to include a timestamp
#' @export
log_warning <- function(message, timestamp = TRUE) {
  if (timestamp) {
    timestamp_str <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    cat(timestamp_str, "WARNING:", message, "\n")
  } else {
    cat("WARNING:", message, "\n")
  }
}

#' Log an error message to the console
#'
#' @param message The error message to log
#' @param timestamp Whether to include a timestamp
#' @param exit Whether to exit the script with an error code
#' @param exit_code The exit code to use if exiting
#' @export
log_error <- function(message, timestamp = TRUE, exit = FALSE, exit_code = 1) {
  if (timestamp) {
    timestamp_str <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    cat(timestamp_str, "ERROR:", message, "\n", file = stderr())
  } else {
    cat("ERROR:", message, "\n", file = stderr())
  }
  
  if (exit) {
    quit(status = exit_code)
  }
}

#' Log details of an error, including the traceback if available
#'
#' @param error The error object
#' @param context Additional context about where the error occurred
#' @param exit Whether to exit the script with an error code
#' @param exit_code The exit code to use if exiting
#' @export
log_error_details <- function(error, context = NULL, exit = FALSE, exit_code = 1) {
  error_message <- conditionMessage(error)
  
  if (!is.null(context)) {
    error_message <- paste0(context, ": ", error_message)
  }
  
  log_error(error_message, timestamp = TRUE, exit = FALSE)
  
  # Log the traceback if available
  if (exists("traceback") && is.function(traceback)) {
    cat("Traceback:\n", file = stderr())
    tryCatch({
      traceback(x = NULL)
    }, error = function(e) {
      cat("Unable to generate traceback\n", file = stderr())
    })
  }
  
  if (exit) {
    quit(status = exit_code)
  }
}

# --- Logging functions ---

# Debug mode control
DEBUG_MODE <- FALSE

# Function to enable debug mode
enable_debug_mode <- function(args) {
  assign("DEBUG_MODE", TRUE, envir = .GlobalEnv)
  log_message("Debug mode enabled")
  
  if (!is.null(args$verbose) && args$verbose) {
    log_message("Verbose logging enabled")
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
    log_message("Started new log file")
  }
}

# Formatted logging function with log levels
log_message_formatted <- function(message, level = "INFO") {
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
    log_message_formatted(message, "DEBUG")
    
    # If data object is provided, print its summary
    if (!is.null(data)) {
      if (is.data.frame(data)) {
        log_message_formatted(glue::glue("Data dimensions: {nrow(data)} rows x {ncol(data)} columns"), "DEBUG")
        if (nrow(data) > 0) {
          class_info <- sapply(data, class)
          log_message_formatted("Column types: ", "DEBUG")
          print(class_info)
          
          # Print a few rows for inspection
          log_message_formatted("Preview of data:", "DEBUG")
          print(head(data, 5))
        }
      } else if (is.matrix(data)) {
        log_message_formatted(glue::glue("Matrix dimensions: {nrow(data)} rows x {ncol(data)} columns"), "DEBUG")
        log_message_formatted("Matrix preview:", "DEBUG")
        print(head(data, 5))
      } else if (is.list(data)) {
        log_message_formatted(glue::glue("List with {length(data)} elements"), "DEBUG")
        log_message_formatted(glue::glue("List names: {paste(names(data), collapse=', ')}"), "DEBUG")
      } else {
        log_message_formatted("Object summary:", "DEBUG")
        print(summary(data))
      }
    }
  }
} 