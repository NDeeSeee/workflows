#!/usr/bin/env Rscript

# --- Utility functions ---

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