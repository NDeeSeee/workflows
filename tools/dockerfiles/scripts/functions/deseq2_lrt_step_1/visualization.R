#!/usr/bin/env Rscript

# --- Visualization functions ---

# Generate PCA plot
generate_pca_plot <- function(dds, metadata_df, args) {
  log_message("Generating PCA plot...", "INFO")
  
  # Perform variance stabilizing transformation
  vsd <- DESeq2::vst(dds, blind = FALSE)
  
  # Extract PCA data
  pca_data <- DESeq2::plotPCA(vsd, intgroup = colnames(metadata_df), returnData = TRUE)
  
  # Calculate percent variance for axis labels
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  # Initialize base plot
  pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2)) +
    ggplot2::labs(
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      title = "Principal Component Analysis",
      subtitle = "Based on variance stabilized expression data"
    )
  
  # Determine colors and shapes based on metadata
  if (ncol(metadata_df) > 0) {
    # First categorical column for color
    color_col <- names(which(sapply(metadata_df, function(x) is.factor(x) || is.character(x))))[1]
    
    if (!is.null(color_col)) {
      pca_plot <- pca_plot + 
        ggplot2::aes(color = !!ggplot2::sym(color_col)) +
        ggplot2::labs(color = color_col)
    }
    
    # Second categorical column for shape (if available)
    shape_cols <- names(which(sapply(metadata_df, function(x) is.factor(x) || is.character(x))))
    
    if (length(shape_cols) > 1) {
      shape_col <- shape_cols[2]
      pca_plot <- pca_plot + 
        ggplot2::aes(shape = !!ggplot2::sym(shape_col)) +
        ggplot2::labs(shape = shape_col)
    }
  }
  
  # Add points and labels
  pca_plot <- pca_plot +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = name), vjust = 1.5, size = 3) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
      legend.position = "right"
    )
  
  # If batch correction was applied, add note to the plot
  if (!is.null(args$batchcorrection) && args$batchcorrection != "none") {
    pca_plot <- pca_plot +
      ggplot2::labs(caption = paste0("Batch correction method: ", args$batchcorrection))
  }
  
  return(pca_plot)
}

# Generate MA plot
generate_ma_plot <- function(lrt_results) {
  log_message("Generating MA plot...", "INFO")
  
  # Convert DESeq2 results to data frame
  lrt_df <- as.data.frame(lrt_results)
  
  # Add column for significance
  lrt_df$significant <- ifelse(!is.na(lrt_df$padj) & lrt_df$padj < 0.1, "Significant", "Not Significant")
  
  # For LRT analysis, log2FoldChange might not be directly available in the same way
  # Adjust the plot accordingly
  if ("log2FoldChange" %in% colnames(lrt_df)) {
    ma_plot <- ggplot2::ggplot(lrt_df, ggplot2::aes(x = baseMean, y = log2FoldChange, color = significant)) +
      ggplot2::geom_point(alpha = 0.6, size = 1) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(
        x = "Mean Expression",
        y = "Log2 Fold Change",
        title = "MA Plot",
        subtitle = "Log2 fold changes vs mean expression",
        color = "FDR < 0.1"
      ) +
      ggplot2::scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray50")) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic")
      )
  } else {
    # If log2FoldChange isn't available, create a basic plot showing significant genes by mean expression
    ma_plot <- ggplot2::ggplot(lrt_df, ggplot2::aes(x = baseMean, y = stat, color = significant)) +
      ggplot2::geom_point(alpha = 0.6, size = 1) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(
        x = "Mean Expression",
        y = "LRT Statistic",
        title = "Expression Plot",
        subtitle = "Test statistic vs mean expression",
        color = "FDR < 0.1"
      ) +
      ggplot2::scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray50")) +
      ggplot2::theme_classic() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic")
      )
  }
  
  return(ma_plot)
}

# Generate p-value histogram
generate_pvalue_histogram <- function(lrt_results) {
  log_message("Generating p-value histogram...", "INFO")
  
  # Convert to data frame
  lrt_df <- as.data.frame(lrt_results)
  
  # Create combined plot (p-value and adjusted p-value)
  p_values <- data.frame(
    value = c(lrt_df$pvalue, lrt_df$padj),
    type = c(rep("P-value", length(lrt_df$pvalue)), rep("Adjusted P-value", length(lrt_df$padj)))
  )
  
  # Filter out NAs
  p_values <- p_values[!is.na(p_values$value), ]
  
  # Create histogram
  p_hist <- ggplot2::ggplot(p_values, ggplot2::aes(x = value, fill = type)) +
    ggplot2::geom_histogram(bins = 30, position = "dodge", alpha = 0.7) +
    ggplot2::labs(
      x = "P-value",
      y = "Count",
      title = "P-value Distribution",
      subtitle = "Distribution of raw and adjusted p-values",
      fill = "Type"
    ) +
    ggplot2::scale_fill_manual(values = c("P-value" = "blue", "Adjusted P-value" = "red")) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic")
    )
  
  return(p_hist)
}

# Generate heatmap of top genes
generate_heatmap <- function(lrt_results, norm_counts, metadata_df, output_prefix, args) {
  log_message("Generating heatmap of top genes...", "INFO")
  
  # Get top genes by adjusted p-value
  top_n <- min(50, nrow(lrt_results))  # Cap at 50 genes or fewer if less are available
  
  top_genes <- rownames(lrt_results)[order(lrt_results$padj)][1:top_n]
  top_genes <- top_genes[!is.na(top_genes)]  # Remove NA values
  
  # Skip if no significant genes
  if (length(top_genes) == 0) {
    log_message("No significant genes to create heatmap", "WARNING")
    return(NULL)
  }
  
  # Extract counts for top genes
  top_counts <- norm_counts[top_genes, ]
  
  # Apply scaling (z-score normalization by row)
  scaled_counts <- t(scale(t(top_counts)))
  
  # Create annotation data for columns based on metadata
  col_annotation <- as.data.frame(metadata_df)
  
  # Generate heatmap with annotations
  pheatmap::pheatmap(
    scaled_counts,
    annotation_col = col_annotation,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 8,
    main = "Heatmap of Top Genes by Adjusted P-value",
    filename = paste0(output_prefix, "_top_genes_heatmap.png"),
    width = 10,
    height = 12
  )
  
  log_message(glue::glue("Generated heatmap with {length(top_genes)} top genes"), "SUCCESS")
} 