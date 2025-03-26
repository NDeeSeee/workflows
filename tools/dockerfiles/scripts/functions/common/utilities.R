#!/usr/bin/env Rscript

# --- Common Utility Functions ---
# This file contains utility functions used across multiple DESeq analysis scripts.
# These functions should be generic and not specific to any particular analysis workflow.

#' Track memory usage during script execution
#' 
#' @param label Label for the memory usage report
#' @return None, prints memory usage message
report_memory_usage <- function(label = "") {
  gc(verbose = FALSE)
  
  if (requireNamespace("pryr", quietly = TRUE)) {
    mem_used <- pryr::mem_used()
    message(paste0("[Memory] ", label, ": ", round(mem_used / 1024^2, 1), " MB"))
  } else {
    message(paste0("[Memory] ", label, ": Not available (pryr package not installed)"))
  }
}

#' Source a file with fallback paths
#' 
#' @param filepath Relative path to file
#' @param absolute_path Absolute path to file (fallback)
#' @return Result of source() call
source_with_fallback <- function(filepath, absolute_path = NULL) {
  # Try absolute path first if provided
  if (!is.null(absolute_path) && file.exists(absolute_path)) {
    message(paste("Sourcing from absolute path:", absolute_path))
    return(source(absolute_path))
  }
  
  # Try relative path
  if (file.exists(filepath)) {
    message(paste("Sourcing from relative path:", filepath))
    return(source(filepath))
  }
  
  # Try a standard Docker path
  docker_path <- file.path("/usr/local/bin", filepath)
  if (file.exists(docker_path)) {
    message(paste("Sourcing from Docker path:", docker_path))
    return(source(docker_path))
  }
  
  # If all fails, error
  stop(paste("Could not find file to source:", filepath))
}

#' Configure standard plotting theme
#' 
#' @return None, sets global ggplot2 theme
configure_plot_theme <- function() {
  # Set default ggplot2 theme if ggplot2 is loaded
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::theme_set(ggplot2::theme_bw() + 
                      ggplot2::theme(text = ggplot2::element_text(size = 12),
                            axis.text = ggplot2::element_text(size = 10),
                            plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                            legend.title = ggplot2::element_text(size = 12),
                            legend.text = ggplot2::element_text(size = 10),
                            strip.background = ggplot2::element_rect(fill = "lightgray"),
                            strip.text = ggplot2::element_text(face = "bold")))
  }
}

#' Write matrix to GCT file
#' 
#' @param data Matrix of expression data
#' @param file_path Path to the output GCT file
#' @export
write_gct <- function(data, file_path) {
  # Check if data is a matrix or data frame
  if (is.data.frame(data)) {
    data_matrix <- as.matrix(data)
  } else {
    data_matrix <- data
  }
  
  # Extract dimensions
  n_rows <- nrow(data_matrix)
  n_cols <- ncol(data_matrix)
  
  # Ensure row and column names exist
  if (is.null(rownames(data_matrix))) {
    rownames(data_matrix) <- paste0("GENE_", 1:n_rows)
  }
  
  if (is.null(colnames(data_matrix))) {
    colnames(data_matrix) <- paste0("SAMPLE_", 1:n_cols)
  }
  
  # Create GCT file format
  tryCatch({
    # Open the file for writing
    con <- file(file_path, "w")
    
    # Write header
    writeLines("#1.2", con)
    writeLines(paste(n_rows, n_cols, sep = "\t"), con)
    
    # Write column headers
    writeLines(paste(c("NAME", "Description", colnames(data_matrix)), collapse = "\t"), con)
    
    # Write data
    for (i in 1:n_rows) {
      # Replace NA with "na" for GCT format
      row_data <- data_matrix[i, ]
      row_data[is.na(row_data)] <- "na"
      
      # Convert to character and format
      row_data_str <- paste(as.character(row_data), collapse = "\t")
      
      # Write row with name and description
      writeLines(paste(rownames(data_matrix)[i], "na", row_data_str, sep = "\t"), con)
    }
    
    # Close file
    close(con)
    
    message(paste("GCT file written successfully to:", file_path))
    
  }, error = function(e) {
    message(paste("Error writing GCT file:", e$message))
    
    # Ensure file is closed
    if (exists("con") && isOpen(con)) {
      close(con)
    }
    
    # Rethrow error
    stop(e)
  })
}

#' Format p-values with scientific notation
#' 
#' @param pvals Vector of p-values
#' @param sig_digits Number of significant digits
#' @return Formatted p-values as character vector
format_pval <- function(pvals, sig_digits = 3) {
  return(
    ifelse(
      is.na(pvals),
      "NA",
      ifelse(
        pvals < 0.001,
        sprintf("%.1e", pvals),
        sprintf(paste0("%.", sig_digits, "f"), pvals)
      )
    )
  )
}

#' Safe version of log2 that handles zeros and negatives
#' 
#' @param x Numeric vector
#' @return log2 transformed values
safe_log2 <- function(x) {
  # Set small values or zeros to a very small number
  x[x <= 0] <- 1e-10
  return(log2(x))
}

#' Safely create a directory
#' 
#' @param path Directory path to create
#' @return Created path
safe_mkdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  return(path)
}

#' Log messages with timestamps and categories
#' 
#' @param message Message to log
#' @param category Log category/level (INFO, WARNING, ERROR)
#' @return None, prints formatted message
log_message <- function(message, category = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_message <- sprintf("[%s] [%s] %s", timestamp, category, message)
  message(formatted_message)
}

#' Report an error with details and recommendations
#' 
#' @param message Error message
#' @param details Optional error details
#' @param recommendations Optional recommendations to resolve error
#' @return None, reports error and stops execution
report_error <- function(message, details = NULL, recommendations = NULL) {
  log_message(paste("ERROR:", message), "ERROR")
  
  if (!is.null(details)) {
    log_message(paste("Details:", details), "ERROR")
  }
  
  if (!is.null(recommendations)) {
    log_message("Recommendations:", "ERROR")
    if (is.character(recommendations) && length(recommendations) > 0) {
      for (rec in recommendations) {
        log_message(paste(" -", rec), "ERROR")
      }
    }
  }
  
  # Stop execution
  stop(message, call. = FALSE)
}

#' Check file delimiter
#' 
#' Determines whether a file uses comma or tab as a delimiter by analyzing the first few lines.
#' 
#' @param file_path Path to the file to check
#' @return Character string indicating the delimiter ("\t" or ",")
check_file_delimiter <- function(file_path) {
  # Read first few lines to determine delimiter
  first_lines <- readLines(file_path, n = 5)
  
  # Count occurrences of common delimiters
  comma_count <- sum(stringr::str_count(first_lines, ","))
  tab_count <- sum(stringr::str_count(first_lines, "\t"))
  
  # Return appropriate delimiter
  if (tab_count > comma_count) {
    return("\t")
  } else {
    return(",")
  }
}

#' Get file type separator
#' 
#' Determines the separator character based on file extension.
#' 
#' @param file_path Path to the file
#' @return Character string containing the separator character ("," or "\t")
get_file_separator <- function(file_path) {
  # Determine separator based on file extension
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    return(",")
  } else if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
    return("\t")
  } else {
    # Auto-detect for other file types
    return(check_file_delimiter(file_path))
  }
}

#' Error handling wrapper
#' 
#' Function to handle errors in a standardized way
#' 
#' @param expr Expression to evaluate
#' @return Result of expression if successful, NULL otherwise
with_error_handling <- function(expr) {
  tryCatch({
    expr
  }, error = function(e) {
    message(paste("ERROR:", e$message))
    message("Stack trace:")
    message(paste(capture.output(print(sys.calls())), collapse = "\n"))
    NULL
  }, warning = function(w) {
    message(paste("WARNING:", w$message))
    NULL
  })
}

#' Load and validate required libraries
#' 
#' @param required_pkgs Vector of package names
#' @param stop_on_missing Whether to stop execution if packages are missing
#' @return Logical indicating whether all packages were loaded
load_required_packages <- function(required_pkgs, stop_on_missing = TRUE) {
  missing_pkgs <- c()
  
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    } else {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
  
  if (length(missing_pkgs) > 0) {
    err_msg <- paste("Missing required packages:", paste(missing_pkgs, collapse = ", "))
    if (stop_on_missing) {
      stop(err_msg)
    } else {
      warning(err_msg)
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#' Remove file extension from filename
#' 
#' @param filename Filename with extension
#' @return Filename without extension
remove_extension <- function(filename) {
  sub("\\.[^.]*$", "", filename)
}

#' Fix column names to make them compatible with R
#' 
#' @param x Vector of column names
#' @return Vector of cleaned column names
fix_colnames <- function(x) {
  # Replace spaces and special characters
  x <- gsub("[^a-zA-Z0-9_]", "_", x)
  
  # Ensure it starts with a letter or underscore
  x <- gsub("^([0-9])", "X\\1", x)
  
  # Make unique if duplicates exist
  if (any(duplicated(x))) {
    message("Warning: Duplicate column names detected and made unique")
    x <- make.unique(x)
  }
  
  return(x)
}

#' Validate metadata
#' 
#' Checks that metadata contains required columns and correct data types
#' 
#' @param metadata_df Data frame containing sample metadata
#' @param batchcorrection Batch correction method to use
#' @param design_formula Design formula for DESeq2
#' @return Validated metadata data frame
#' @export
validate_metadata <- function(metadata_df, batchcorrection = "none", design_formula = NULL) {
  # Check if metadata has any rows
  if (nrow(metadata_df) == 0) {
    stop("Metadata has no rows")
  }
  
  # Check for batch correction requirements
  if (batchcorrection != "none") {
    if (!"batch" %in% colnames(metadata_df)) {
      warning("Batch correction requested but 'batch' column not found in metadata. Batch correction will be disabled.")
    } else {
      # Convert batch to numeric if it's not already
      if (!is.numeric(metadata_df$batch)) {
        warning("Converting 'batch' column to numeric")
        metadata_df$batch <- as.numeric(as.factor(metadata_df$batch))
      }
    }
  }
  
  # If design formula is provided, validate factors referenced in it
  if (!is.null(design_formula)) {
    formula_vars <- all.vars(design_formula)
    missing_vars <- formula_vars[!formula_vars %in% colnames(metadata_df)]
    
    if (length(missing_vars) > 0) {
      stop(paste("Design formula variables not found in metadata:", paste(missing_vars, collapse=", ")))
    }
    
    # Convert character columns used in design to factors
    for (var in formula_vars) {
      if (is.character(metadata_df[[var]])) {
        metadata_df[[var]] <- as.factor(metadata_df[[var]])
        message(paste("Converted", var, "to factor"))
      }
    }
  }
  
  # Return the validated and possibly modified metadata
  return(metadata_df)
}

#' Validate sample consistency
#' 
#' Checks that sample names in metadata and count data match
#' 
#' @param metadata_df Data frame containing sample metadata
#' @param counts_df Data frame containing count data
#' @return Nothing, stops execution if inconsistencies are found
#' @export
validate_sample_consistency <- function(metadata_df, counts_df) {
  # Get sample names from both dataframes
  metadata_samples <- rownames(metadata_df)
  count_samples <- colnames(counts_df)
  
  # Check if all metadata samples are in count data
  missing_in_counts <- base::setdiff(metadata_samples, count_samples)
  if (length(missing_in_counts) > 0) {
    stop(paste("Samples in metadata but not in count data:", paste(missing_in_counts, collapse=", ")))
  }
  
  # Check if all count data samples are in metadata
  missing_in_metadata <- base::setdiff(count_samples, metadata_samples)
  if (length(missing_in_metadata) > 0) {
    stop(paste("Samples in count data but not in metadata:", paste(missing_in_metadata, collapse=", ")))
  }
  
  # If we got this far, samples are consistent
  message(paste("Sample consistency check passed for", length(metadata_samples), "samples"))
}

#' Format and print all arguments for logging
#' 
#' @param args List of arguments to format
#' @return Formatted string with all arguments
#' @export
print_all_args <- function(args) {
  # Convert arguments to a character vector for printing
  args_str <- capture.output(print(args))
  
  # Format as a string with newlines
  return(paste(args_str, collapse = "\n"))
}

#' Set logger to verbose level
#' 
#' Configures the logger package to use verbose logging level
#' 
#' @return Nothing, modifies global logger settings
#' @export
set_log_level_verbose <- function() {
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_threshold(logger::DEBUG)
    message("Logger configured for verbose output")
  } else {
    message("Warning: logger package not available, using base message() function")
  }
}

#' Log debug level message
#' 
#' @param ... Arguments passed to logger::log_debug or message
#' @export
log_debug <- function(...) {
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_debug(...)
  } else {
    message("[DEBUG] ", ...)
  }
}

#' Log info level message
#' 
#' @param ... Arguments passed to logger::log_info or message
#' @export
log_info <- function(...) {
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_info(...)
  } else {
    message("[INFO] ", ...)
  }
}

#' Log error level message
#' 
#' @param ... Arguments passed to logger::log_error or message
#' @export
log_error <- function(...) {
  if (requireNamespace("logger", quietly = TRUE)) {
    logger::log_error(...)
  } else {
    message("[ERROR] ", ...)
  }
} 