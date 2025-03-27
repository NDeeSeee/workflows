#!/usr/bin/env Rscript
#
# Output utility functions for DESeq2 workflows
# Contains minimal functions needed to ensure proper output formats
#

#' Create the exact file name expected by CWL
#'
#' @param prefix User-provided output prefix
#' @param stem Middle part of the filename
#' @param extension File extension (without the dot)
#' @return Full filename
#' @export
get_output_filename <- function(prefix, stem, extension) {
  # Clean prefix (remove trailing slashes)
  prefix <- gsub("/*$", "", prefix)
  
  if (stem == "") {
    return(paste0(prefix, ".", extension))
  } else {
    return(paste0(prefix, "_", stem, ".", extension))
  }
}

#' Save a plot in both PNG and PDF formats
#'
#' @param plot A ggplot2 object or function that creates a plot
#' @param prefix Output prefix
#' @param stem Middle part of the filename
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return List of created filenames
#' @export
save_plot_png_pdf <- function(plot, prefix, stem, width = 8, height = 6) {
  results <- list()
  
  # PNG output
  png_file <- get_output_filename(prefix, stem, "png")
  pdf_file <- get_output_filename(prefix, stem, "pdf")
  
  # For ggplot2 plots
  if (inherits(plot, "ggplot")) {
    # Save PNG
    ggplot2::ggsave(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      dpi = 300,
      units = "in"
    )
    
    # Save PDF
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot,
      width = width,
      height = height,
      device = cairo_pdf
    )
  } else {
    # For base R plots or function calls
    # PNG
    png(png_file, width = width, height = height, units = "in", res = 300)
    if (is.function(plot)) plot() else print(plot)
    dev.off()
    
    # PDF
    pdf(pdf_file, width = width, height = height)
    if (is.function(plot)) plot() else print(plot)
    dev.off()
  }
  
  results$png <- png_file
  results$pdf <- pdf_file
  
  return(results)
}

#' Write a GCT file for GSEA
#'
#' @param data Expression data matrix or data frame (genes in rows, samples in columns)
#' @param file_path Output file path
#' @return The path to the created file
#' @export
write_gct_file <- function(data, file_path) {
  # Make sure it's a matrix or data frame
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input data must be a matrix or data frame")
  }
  
  # Extract dimensions
  n_rows <- nrow(data)
  n_cols <- ncol(data)
  
  # Open connection to file
  con <- file(file_path, "w")
  
  # Write GCT header
  writeLines("#1.2", con)
  writeLines(paste(n_rows, n_cols, sep = "\t"), con)
  
  # Write column headers
  col_headers <- c("NAME", "Description", colnames(data))
  writeLines(paste(col_headers, collapse = "\t"), con)
  
  # Write data
  for (i in 1:n_rows) {
    # Create description (use "na" as a placeholder)
    description <- "na"
    
    # Get data for the row
    row_data <- data[i, ]
    row_name <- rownames(data)[i]
    
    # Format the output line
    row_str <- paste(c(row_name, description, as.character(row_data)), collapse = "\t")
    writeLines(row_str, con)
  }
  
  close(con)
  message(paste("Created GCT file:", file_path))
  
  return(file_path)
}

#' Create a CLS file for GSEA
#'
#' @param sample_classes Vector of sample class labels (must match order of columns in expression data)
#' @param file_path Output file path
#' @return The path to the created file
#' @export
write_cls_file <- function(sample_classes, file_path) {
  # Get unique classes
  unique_classes <- unique(sample_classes)
  n_classes <- length(unique_classes)
  n_samples <- length(sample_classes)
  
  # Create header line
  header <- paste(n_samples, n_classes, 1, sep = " ")
  
  # Create class names line
  class_names <- paste("#", paste(unique_classes, collapse = " "))
  
  # Create sample assignments line
  # Convert each class to its index (0-based)
  class_indices <- sapply(sample_classes, function(x) {
    which(unique_classes == x) - 1
  })
  
  assignments <- paste(class_indices, collapse = " ")
  
  # Write to file
  writeLines(c(header, class_names, assignments), file_path)
  message(paste("Created CLS file:", file_path))
  
  return(file_path)
}

#' Create a simple markdown summary file
#'
#' @param content Character vector of content lines
#' @param file_path Output file path
#' @return The path to the created file
#' @export
write_markdown_summary <- function(content, file_path) {
  writeLines(content, file_path)
  message(paste("Created markdown summary:", file_path))
  
  return(file_path)
}

#' Check if all expected outputs exist
#'
#' @param expected_files List of filenames that should exist
#' @param stop_on_missing Whether to stop execution if files are missing
#' @return Logical indicating if all files exist
#' @export
check_outputs <- function(expected_files, stop_on_missing = TRUE) {
  missing_files <- character(0)
  
  for (file in expected_files) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) > 0) {
    msg <- paste0("Missing expected output files: ", 
                 paste(missing_files, collapse = ", "))
    
    if (stop_on_missing) {
      stop(msg)
    } else {
      warning(msg)
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#' Get list of expected outputs for a specific workflow
#'
#' @param workflow_name Name of the workflow: "deseq-advanced", "lrt-step-1", or "lrt-step-2"
#' @param prefix Output prefix
#' @return Character vector of expected output filenames
#' @export
get_expected_outputs <- function(workflow_name, prefix) {
  # Clean prefix
  prefix <- gsub("/*$", "", prefix)
  
  if (workflow_name == "deseq-advanced") {
    return(c(
      paste0(prefix, "_report.tsv"),
      paste0(prefix, "_summary.md"),
      paste0(prefix, "_counts_all.gct"),
      paste0(prefix, "_counts_filtered.gct"),
      paste0(prefix, "_phenotypes.cls"),
      paste0(prefix, "_ma_plot.png"),
      paste0(prefix, "_expression_heatmap.png"),
      paste0(prefix, "_pca_plot.png"),
      paste0(prefix, "_ma_plot.pdf"),
      paste0(prefix, "_expression_heatmap.pdf"),
      paste0(prefix, "_pca_plot.pdf"),
      paste0(prefix, "_mds_plot.html")
    ))
  } else if (workflow_name == "lrt-step-1") {
    return(c(
      paste0(prefix, "_contrasts_table.tsv"),
      paste0(prefix, "_contrasts.rds"),
      paste0(prefix, "_gene_exp_table.tsv"),
      paste0(prefix, "_mds_plot.html"),
      paste0(prefix, "_mds_plot_corrected.html"),
      paste0(prefix, "_counts_all.gct"),
      paste0(prefix, "_counts_filtered.gct"),
      paste0(prefix, "_lrt_result.md"),
      "alignment_stats_barchart.png"
    ))
  } else if (workflow_name == "lrt-step-2") {
    return(c(
      paste0(prefix, "_gene_exp_table.tsv"),
      "mds_plot.html",
      "counts_all.gct",
      "counts_filtered.gct"
    ))
  } else {
    stop(paste("Unknown workflow:", workflow_name))
  }
} 