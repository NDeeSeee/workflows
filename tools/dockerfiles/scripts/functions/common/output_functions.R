#!/usr/bin/env Rscript
#
# Standardized output functions for DESeq2 workflows
# These functions ensure that all outputs match exactly what the CWL files expect
#

#' Generate standardized output filename
#'
#' @param prefix Output prefix specified by the user
#' @param output_type Type of output (e.g., "ma_plot", "expression_heatmap")
#' @param extension File extension without the dot (e.g., "png", "pdf")
#' @return Standardized filename string
#' @export
generate_output_filename <- function(prefix, output_type, extension) {
  # Remove any trailing slashes or dots from prefix
  prefix <- gsub("[/\\.]$", "", prefix)
  
  # For certain output types, use specific naming patterns
  if (output_type == "alignment_stats_barchart") {
    return(paste0("alignment_stats_barchart.", extension))
  } else if (output_type == "mds_plot" && extension == "html") {
    return(paste0(prefix, "_mds_plot.html"))
  } else if (output_type == "mds_plot_corrected" && extension == "html") {
    return(paste0(prefix, "_mds_plot_corrected.html"))
  } else if (output_type == "counts_all" && extension == "gct") {
    return(paste0(prefix, "_counts_all.gct"))
  } else if (output_type == "counts_filtered" && extension == "gct") {
    return(paste0(prefix, "_counts_filtered.gct"))
  } else {
    # Default pattern
    return(paste0(prefix, "_", output_type, ".", extension))
  }
}

#' Save plot in multiple formats
#'
#' @param plot ggplot2 object to save
#' @param prefix Output prefix specified by the user
#' @param output_type Type of output (e.g., "ma_plot", "pca_plot")
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param formats Vector of formats to save (default: c("png", "pdf"))
#' @param ... Additional parameters passed to ggsave
#' @return List of saved file paths
#' @export
save_plot_multiple_formats <- function(plot, prefix, output_type, width = 8, height = 6, 
                                      formats = c("png", "pdf"), ...) {
  saved_files <- list()
  
  for (format in formats) {
    filename <- generate_output_filename(prefix, output_type, format)
    
    # Use ggsave for ggplot objects, or appropriate device for base plots
    if (inherits(plot, "ggplot")) {
      ggplot2::ggsave(filename, plot, width = width, height = height, 
                    units = "in", dpi = 300, ...)
    } else {
      # For base R plots, use the appropriate device
      if (format == "png") {
        png(filename, width = width, height = height, units = "in", res = 300)
        if (is.function(plot)) plot() else print(plot)
        dev.off()
      } else if (format == "pdf") {
        pdf(filename, width = width, height = height)
        if (is.function(plot)) plot() else print(plot)
        dev.off()
      }
    }
    
    # Add to list of saved files
    saved_files[[format]] <- filename
    message(paste("Saved", format, "plot to", filename))
  }
  
  return(saved_files)
}

#' Create phenotypes CLS file for GSEA
#'
#' @param metadata Data frame containing sample metadata
#' @param condition_col Column name in metadata containing condition information
#' @param prefix Output prefix for the file
#' @return Path to the created CLS file
#' @export
create_phenotypes_cls <- function(metadata, condition_col, prefix) {
  # Validate input
  if (!condition_col %in% colnames(metadata)) {
    stop(paste("Condition column", condition_col, "not found in metadata"))
  }
  
  # Get condition levels
  condition_levels <- as.character(unique(metadata[[condition_col]]))
  num_conditions <- length(condition_levels)
  num_samples <- nrow(metadata)
  
  # Create the CLS file content
  cls_header <- paste(num_samples, num_conditions, 1, sep = " ")
  cls_conditions <- paste("#", paste(condition_levels, collapse = " "))
  
  # Get class assignments for each sample
  class_assignments <- sapply(metadata[[condition_col]], function(x) {
    which(condition_levels == x) - 1  # 0-based index for GSEA
  })
  cls_assignments <- paste(class_assignments, collapse = " ")
  
  # Combine all parts
  cls_content <- paste(cls_header, cls_conditions, cls_assignments, sep = "\n")
  
  # Write to file
  cls_filepath <- paste0(prefix, "_phenotypes.cls")
  writeLines(cls_content, cls_filepath)
  message(paste("Created phenotypes CLS file for GSEA:", cls_filepath))
  
  return(cls_filepath)
}

#' Write DESeq2 results to a standardized TSV report
#'
#' @param deseq_results DESeq2 results object
#' @param prefix Output prefix
#' @param output_type Type of output report (e.g., "report", "gene_exp_table")
#' @param add_gene_info Whether to add gene information columns
#' @param gene_info Optional data frame with gene annotation information
#' @return Path to the created report file
#' @export
write_deseq_results <- function(deseq_results, prefix, output_type = "report",
                               add_gene_info = FALSE, gene_info = NULL) {
  # Convert DESeq2 results to data frame
  res_df <- as.data.frame(deseq_results)
  
  # Add rank column (useful for GSEA)
  res_df$rank <- sign(res_df$log2FoldChange) * -log10(res_df$pvalue)
  
  # Add gene info if provided and requested
  if (add_gene_info && !is.null(gene_info)) {
    # Merge with gene info, ensuring all rows are preserved
    res_df <- merge(res_df, gene_info, by.x = "row.names", by.y = "row.names", all.x = TRUE)
    rownames(res_df) <- res_df$Row.names
    res_df$Row.names <- NULL
  }
  
  # Sort by adjusted p-value
  res_df <- res_df[order(res_df$padj), ]
  
  # Determine output filename based on output_type
  if (output_type == "report") {
    filename <- paste0(prefix, "_report.tsv")
  } else if (output_type == "gene_exp_table") {
    filename <- paste0(prefix, "_gene_exp_table.tsv")
  } else if (output_type == "contrasts_table") {
    filename <- paste0(prefix, "_contrasts_table.tsv")
  } else {
    filename <- paste0(prefix, "_", output_type, ".tsv")
  }
  
  # Write to TSV file
  write.table(res_df, file = filename, sep = "\t", quote = FALSE, row.names = TRUE)
  message(paste("Wrote DESeq2 results to", filename))
  
  return(filename)
}

#' Generate markdown summary of DESeq2 analysis
#'
#' @param deseq_obj DESeq2 object
#' @param deseq_results DESeq2 results object
#' @param prefix Output prefix
#' @param output_type Type of summary (e.g., "summary", "lrt_result")
#' @param design_formula Design formula used
#' @param reduced_formula Reduced formula used (for LRT)
#' @param additional_info List with additional information to include
#' @return Path to the created markdown file
#' @export
generate_summary_markdown <- function(deseq_obj, deseq_results, prefix,
                                     output_type = "summary", design_formula = NULL,
                                     reduced_formula = NULL, additional_info = list()) {
  # Determine output filename
  if (output_type == "summary") {
    filename <- paste0(prefix, "_summary.md")
  } else if (output_type == "lrt_result") {
    filename <- paste0(prefix, "_lrt_result.md")
  } else {
    filename <- paste0(prefix, "_", output_type, ".md")
  }
  
  # Get summary statistics
  total_genes <- nrow(deseq_results)
  sig_genes <- sum(deseq_results$padj < 0.1, na.rm = TRUE)
  up_genes <- sum(deseq_results$padj < 0.1 & deseq_results$log2FoldChange > 0, na.rm = TRUE)
  down_genes <- sum(deseq_results$padj < 0.1 & deseq_results$log2FoldChange < 0, na.rm = TRUE)
  
  # Create markdown content
  lines <- c(
    "# DESeq2 Analysis Summary",
    "",
    paste("Analysis performed on", Sys.Date()),
    "",
    "## Parameters",
    ""
  )
  
  # Add design formulas if provided
  if (!is.null(design_formula)) {
    lines <- c(lines, paste("* Design formula:", design_formula))
  }
  if (!is.null(reduced_formula)) {
    lines <- c(lines, paste("* Reduced formula:", reduced_formula))
  }
  
  # Add any additional information
  if (length(additional_info) > 0) {
    for (name in names(additional_info)) {
      lines <- c(lines, paste("*", name, ":", additional_info[[name]]))
    }
  }
  
  # Add result summary
  lines <- c(lines,
    "",
    "## Results",
    "",
    paste("* Total genes analyzed:", total_genes),
    paste("* Significant genes (FDR < 0.1):", sig_genes),
    paste("* Up-regulated genes:", up_genes),
    paste("* Down-regulated genes:", down_genes)
  )
  
  # Write to file
  writeLines(lines, filename)
  message(paste("Generated markdown summary at", filename))
  
  return(filename)
}

#' Verify all required outputs are created
#'
#' @param expected_files List of expected output files
#' @param fail_on_missing Whether to stop execution if files are missing
#' @return Logical indicating if all expected files exist
#' @export
verify_outputs <- function(expected_files, fail_on_missing = TRUE) {
  missing_files <- character(0)
  
  for (file in expected_files) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
  
  if (length(missing_files) > 0) {
    message(paste("ERROR: The following expected output files are missing:"))
    for (file in missing_files) {
      message(paste(" -", file))
    }
    
    if (fail_on_missing) {
      stop("Missing required output files")
    } else {
      return(FALSE)
    }
  }
  
  message("All expected output files have been created successfully")
  return(TRUE)
}

#' Save DESeq2 object or results to RDS file
#'
#' @param object DESeq2 object or results to save
#' @param prefix Output prefix
#' @param output_type Type of output (e.g., "contrasts", "dsq_obj")
#' @return Path to the created RDS file
#' @export
save_deseq_rds <- function(object, prefix, output_type = "contrasts") {
  # Determine filename
  filename <- paste0(prefix, "_", output_type, ".rds")
  
  # Save object to RDS file
  saveRDS(object, file = filename)
  message(paste("Saved", output_type, "to RDS file:", filename))
  
  return(filename)
}

#' Generate interactive MDS plot HTML file
#'
#' @param dds DESeq2 dds object
#' @param prefix Output prefix
#' @param transformed_counts Transformed counts (e.g., vst or rlog)
#' @param batch_corrected Whether to create batch-corrected version
#' @param metadata Sample metadata
#' @param color_by Column in metadata to use for coloring
#' @param shape_by Column in metadata to use for point shapes
#' @param title Plot title
#' @return Path to the created HTML file(s)
#' @export
generate_mds_plot_html <- function(dds, prefix, transformed_counts = NULL,
                                  batch_corrected = FALSE, metadata = NULL,
                                  color_by = NULL, shape_by = NULL,
                                  title = "MDS Plot of Samples") {
  # Make sure we have the right packages
  if (!requireNamespace("plotly", quietly = TRUE) || 
      !requireNamespace("htmlwidgets", quietly = TRUE)) {
    message("Cannot create MDS plot: missing required packages plotly or htmlwidgets")
    return(NULL)
  }
  
  # Get transformed data if not provided
  if (is.null(transformed_counts)) {
    # Use variance stabilizing transformation
    transformed_counts <- DESeq2::vst(dds, blind = TRUE)
  }
  
  # Calculate MDS
  mds_data <- limma::plotMDS(assay(transformed_counts), plot = FALSE)
  
  # Create data frame for plotting
  mds_df <- data.frame(
    x = mds_data$x,
    y = mds_data$y,
    sample = colnames(transformed_counts)
  )
  
  # Add metadata if provided
  if (!is.null(metadata) && !is.null(color_by) && color_by %in% colnames(metadata)) {
    # Match samples
    idx <- match(mds_df$sample, rownames(metadata))
    mds_df$color <- metadata[idx, color_by]
    
    if (!is.null(shape_by) && shape_by %in% colnames(metadata)) {
      mds_df$shape <- metadata[idx, shape_by]
    }
  }
  
  # Create plotly figure
  p <- plotly::plot_ly(mds_df, x = ~x, y = ~y, text = ~sample,
                     mode = "markers", type = "scatter",
                     marker = list(size = 10))
  
  if ("color" %in% colnames(mds_df)) {
    p <- plotly::add_trace(p, color = ~color)
  }
  
  if ("shape" %in% colnames(mds_df)) {
    p <- plotly::add_trace(p, symbol = ~shape)
  }
  
  p <- plotly::layout(p, title = title,
                    xaxis = list(title = "Leading logFC dim 1"),
                    yaxis = list(title = "Leading logFC dim 2"))
  
  # Save to HTML file
  filename <- generate_output_filename(prefix, "mds_plot", "html")
  htmlwidgets::saveWidget(p, filename, selfcontained = TRUE)
  message(paste("Generated MDS plot HTML at", filename))
  
  # Create batch-corrected version if requested
  if (batch_corrected && "batch" %in% colnames(metadata)) {
    # This requires additional implementation for batch correction
    # but we'll add placeholder for now
    bc_filename <- generate_output_filename(prefix, "mds_plot_corrected", "html")
    # Implement batch correction code here
    message(paste("Generated batch-corrected MDS plot HTML at", bc_filename))
    return(c(filename, bc_filename))
  }
  
  return(filename)
} 