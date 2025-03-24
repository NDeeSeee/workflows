#!/usr/bin/env Rscript

# --- Utility functions ---

# Function to determine file type based on extension
get_file_type <- function(file_path) {
  # Get file extension
  ext <- tolower(tools::file_ext(file_path))
  
  # Return appropriate delimiter based on extension
  if (ext == "csv") {
    return(",")
  } else if (ext %in% c("tsv", "txt")) {
    return("\t")
  } else {
    # Default to tab if extension is unknown
    log_message(glue::glue("Warning: Unknown file extension '{ext}', defaulting to tab delimiter"), "WARNING")
    return("\t")
  }
}

# Function to check file delimiter
check_file_delimiter <- function(file_path) {
  # Get the first line of the file
  first_line <- readLines(file_path, n = 1)
  
  # Count occurrences of common delimiters
  comma_count <- stringr::str_count(first_line, ",")
  tab_count <- stringr::str_count(first_line, "\t")
  
  # Determine expected delimiter based on file extension
  expected_delimiter <- get_file_type(file_path)
  expected_name <- ifelse(expected_delimiter == ",", "comma", "tab")
  
  # Check if the delimiter usage matches the file extension
  if (expected_delimiter == "," && tab_count > comma_count) {
    report_error(
      "File format mismatch",
      details = c(
        "The file appears to be tab-delimited, but has a .csv extension.",
        paste0("First line contains ", tab_count, " tabs and ", comma_count, " commas.")
      ),
      recommendations = c(
        "Rename the file with a .tsv extension",
        "OR convert the file to use comma delimiters to match the .csv extension"
      )
    )
  } else if (expected_delimiter == "\t" && comma_count > tab_count) {
    report_error(
      "File format mismatch",
      details = c(
        "The file appears to be comma-delimited, but has a .tsv extension.",
        paste0("First line contains ", comma_count, " commas and ", tab_count, " tabs.")
      ),
      recommendations = c(
        "Rename the file with a .csv extension",
        "OR convert the file to use tab delimiters to match the .tsv extension"
      )
    )
  }
  
  # Warn if first row has no delimiters
  if (comma_count == 0 && tab_count == 0) {
    log_message(
      "Warning: First line of file contains no delimiters. This may indicate a formatting problem.",
      "WARNING"
    )
  }
}

# Function to clean sample names
clean_sample_names <- function(sample_names) {
  # Replace spaces and special characters with underscores
  cleaned_names <- gsub("[^a-zA-Z0-9]", "_", sample_names)
  
  # Remove duplicate underscores
  cleaned_names <- gsub("_+", "_", cleaned_names)
  
  # Remove leading/trailing underscores
  cleaned_names <- gsub("^_|_$", "", cleaned_names)
  
  return(cleaned_names)
}

# Utility function to remove row names from a data frame
remove_rownames <- function(df) {
  rownames(df) <- NULL
  return(df)
}

# Configure plot theme for consistency
configure_plot_theme <- function() {
  # Set default theme for all plots
  ggplot2::theme_set(
    ggplot2::theme_classic() + 
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      legend.position = "right"
    )
  )
} 