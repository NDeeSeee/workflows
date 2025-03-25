#!/usr/bin/env Rscript

# --- Utility functions ---

# Function to track memory usage
report_memory_usage <- function(label = "") {
  gc(verbose = FALSE)
  mem_used <- pryr::mem_used()
  message(paste0("[Memory] ", label, ": ", round(mem_used / 1024^2, 1), " MB"))
}

# Helper function to source files with fallback paths
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
  
  # If we get here, try a standard Docker path
  docker_path <- file.path("/usr/local/bin", filepath)
  if (file.exists(docker_path)) {
    message(paste("Sourcing from Docker path:", docker_path))
    return(source(docker_path))
  }
  
  # If all fails, error
  stop(paste("Could not find file to source:", filepath))
}

# Configure standard plotting theme
configure_plot_theme <- function() {
  # Set default ggplot2 theme if ggplot2 is loaded
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::theme_set(ggplot2::theme_bw())
  }
}

# Function to write data in GCT format
write_gct <- function(expression_data, output_file) {
  # Check if cmapR is available, otherwise use our own implementation
  if (requireNamespace("cmapR", quietly = TRUE)) {
    # Create a GCT object
    gct_obj <- cmapR::new.gct(expression_data)
    # Write the GCT file
    cmapR::write.gct(gct_obj, output_file)
  } else {
    # Fallback implementation
    write_gct_fallback(expression_data, output_file)
  }
}

# Fallback function to write GCT files when cmapR is not available
write_gct_fallback <- function(expression_data, output_file) {
  # Make sure we have rownames and colnames
  if (is.null(rownames(expression_data))) {
    rownames(expression_data) <- paste0("gene_", 1:nrow(expression_data))
  }
  if (is.null(colnames(expression_data))) {
    colnames(expression_data) <- paste0("sample_", 1:ncol(expression_data))
  }
  
  # Create row descriptions (using row names as a placeholder)
  row_descriptions <- rownames(expression_data)
  
  # Create header
  header <- c(
    paste("#1.2"),
    paste(nrow(expression_data), ncol(expression_data), sep="\t")
  )
  
  # Create column names line
  col_header <- c("NAME", "Description", colnames(expression_data))
  
  # Open the file for writing
  con <- file(output_file, "w")
  
  # Write header
  writeLines(header, con)
  
  # Write column names
  writeLines(paste(col_header, collapse="\t"), con)
  
  # Write data rows
  for (i in 1:nrow(expression_data)) {
    row_data <- c(
      rownames(expression_data)[i],
      row_descriptions[i],
      expression_data[i, ]
    )
    writeLines(paste(row_data, collapse="\t"), con)
  }
  
  # Close the file
  close(con)
  
  message(paste("GCT file written to:", output_file))
}

# Function to format p-values with scientific notation
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

# Safe version of log2 that handles zeros and negatives
safe_log2 <- function(x) {
  # Set small values or zeros to a very small number
  x[x <= 0] <- 1e-10
  return(log2(x))
}

# Function to safely create a directory
safe_mkdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  return(path)
}

# Function to log messages with timestamps and categories
log_message <- function(message, category = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_message <- sprintf("[%s] [%s] %s", timestamp, category, message)
  message(formatted_message)
}

# Function to report an error with details
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
#' @export
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
#' @export
get_file_type <- function(file_path) {
  # Get file extension
  ext <- tolower(tools::file_ext(file_path))
  
  # Determine separator based on extension
  if (ext == "csv") {
    return(",")
  } else if (ext == "tsv" || ext == "txt") {
    return("\t")
  } else {
    # Default to comma with a warning for unknown extensions
    warning(paste("Unknown file extension:", ext, "defaulting to comma separator"))
    return(",")
  }
}

#' Clean sample names 
#' 
#' Standardizes sample names by removing special characters and spaces
#' 
#' @param sample_names Vector of sample names to clean
#' @return Vector of cleaned sample names
#' @export
clean_sample_names <- function(sample_names) {
  # Replace spaces with underscores
  sample_names <- gsub(" ", "_", sample_names)
  
  # Remove special characters except underscores and alphanumeric
  sample_names <- gsub("[^a-zA-Z0-9_]", "", sample_names)
  
  # Ensure names are unique
  if (any(duplicated(sample_names))) {
    warning("Duplicate sample names found after cleaning. Adding unique suffixes.")
    dupes <- which(duplicated(sample_names))
    for (i in dupes) {
      sample_names[i] <- paste0(sample_names[i], "_", i)
    }
  }
  
  return(sample_names)
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
  missing_in_counts <- setdiff(metadata_samples, count_samples)
  if (length(missing_in_counts) > 0) {
    stop(paste("Samples in metadata but not in count data:", paste(missing_in_counts, collapse=", ")))
  }
  
  # Check if all count data samples are in metadata
  missing_in_metadata <- setdiff(count_samples, metadata_samples)
  if (length(missing_in_metadata) > 0) {
    stop(paste("Samples in count data but not in metadata:", paste(missing_in_metadata, collapse=", ")))
  }
  
  # If we got this far, samples are consistent
  message(paste("Sample consistency check passed for", length(metadata_samples), "samples"))
} 