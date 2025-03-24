#!/usr/bin/env Rscript

# --- Export functions ---

# Export results and visualizations
export_results <- function(deseq_results, expression_df, metadata_df, args, batch_warning, rpkm_filtered_count) {
  log_message("Exporting results and visualizations...", "STEP")
  
  # Create output directory if it doesn't exist
  output_dir <- args$output
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set output prefix
  output_prefix <- file.path(output_dir, args$prefix)
  
  # Save DESeq2 dataset and results as RDS files
  saveRDS(deseq_results$dds, file = paste0(output_prefix, "_dds.rds"))
  saveRDS(deseq_results$lrt_res, file = paste0(output_prefix, "_lrt_results.rds"))
  saveRDS(expression_df, file = paste0(output_prefix, "_expression_data.rds"))
  saveRDS(metadata_df, file = paste0(output_prefix, "_metadata.rds"))
  saveRDS(deseq_results$contrasts, file = paste0(output_prefix, "_contrasts.rds"))
  
  # Save batch correction info
  saveRDS(list(
    applied = !is.null(args$batchcorrection) && args$batchcorrection != "none",
    method = args$batchcorrection,
    warning = batch_warning
  ), file = paste0(output_prefix, "_batch_correction_info.rds"))
  
  # Export data table as TSV
  export_table(
    deseq_results$lrt_res,
    expression_df,
    output_prefix,
    args
  )
  
  # Export contrast information
  export_contrasts(
    deseq_results$contrasts,
    output_prefix
  )
  
  # Generate and export visualizations
  export_visualizations(
    deseq_results$dds,
    deseq_results$lrt_res,
    metadata_df,
    deseq_results$normCounts, 
    output_prefix, 
    args
  )
  
  # Export analysis parameters
  export_parameters(
    args, 
    rpkm_filtered_count,
    output_prefix
  )
  
  log_message("Results exported successfully", "SUCCESS")
}

# Export data table as TSV
export_table <- function(lrt_res, expression_df, output_prefix, args) {
  log_message("Exporting data table...", "INFO")
  
  # Convert DESeq2 results to data frame
  lrt_res_df <- as.data.frame(lrt_res)
  
  # Add gene names as column
  lrt_res_df$GeneId <- rownames(lrt_res_df)
  
  # Reorganize columns
  lrt_res_df <- lrt_res_df %>%
    dplyr::select(GeneId, everything())
  
  # Add negative log10 p-values
  lrt_res_df <- lrt_res_df %>%
    dplyr::mutate(
      `-log10(pvalue)` = -log10(pvalue),
      `-log10(padj)` = -log10(padj)
    )
  
  # Save as TSV
  write.table(
    lrt_res_df,
    file = paste0(output_prefix, "_results.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  log_message(glue::glue("Exported {nrow(lrt_res_df)} rows to {output_prefix}_results.tsv"), "SUCCESS")
}

# Export contrast information
export_contrasts <- function(contrasts, output_prefix) {
  log_message("Exporting contrast information...", "INFO")
  
  # Write contrasts to TSV
  write.table(
    contrasts,
    file = paste0(output_prefix, "_contrasts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  log_message(glue::glue("Exported {nrow(contrasts)} contrasts to {output_prefix}_contrasts.tsv"), "SUCCESS")
}

# Export visualizations
export_visualizations <- function(dds, lrt_results, metadata_df, norm_counts, output_prefix, args) {
  log_message("Generating visualizations...", "STEP")
  
  # Generate PCA plot
  pca_plot <- generate_pca_plot(dds, metadata_df, args)
  ggsave(paste0(output_prefix, "_pca.png"), plot = pca_plot, width = 10, height = 8)
  
  # Generate MA plot
  ma_plot <- generate_ma_plot(lrt_results)
  ggsave(paste0(output_prefix, "_ma_plot.png"), plot = ma_plot, width = 10, height = 8)
  
  # Generate p-value distribution
  pval_plot <- generate_pvalue_histogram(lrt_results)
  ggsave(paste0(output_prefix, "_pvalue_hist.png"), plot = pval_plot, width = 10, height = 6)
  
  # Generate heatmap of top genes
  generate_heatmap(lrt_results, norm_counts, metadata_df, output_prefix, args)
  
  # Generate GCT files for visualization in external tools
  export_gct(norm_counts, lrt_results, output_prefix, args)
  
  log_message("Visualizations exported successfully", "SUCCESS")
}

# Export GCT files
export_gct <- function(norm_counts, lrt_results, output_prefix, args) {
  log_message("Exporting GCT files...", "INFO")
  
  # Create row metadata
  row_metadata <- data.frame(
    id = rownames(norm_counts),
    Description = rownames(norm_counts),
    stringsAsFactors = FALSE
  )
  
  # Add columns from DESeq2 results
  lrt_df <- as.data.frame(lrt_results)
  lrt_df$id <- rownames(lrt_df)
  
  # Merge with row metadata
  row_metadata_merged <- dplyr::full_join(
    row_metadata,
    lrt_df %>% dplyr::select(id, baseMean, log2FoldChange, pvalue, padj),
    by = "id"
  )
  
  # Export complete normalized count matrix
  log_message("Exporting complete normalized count matrix...", "INFO")
  cmapR::write_gct(
    ds = cmapR::new_gct(
      mat = as.matrix(norm_counts),
      rid = row_metadata_merged$id,
      rdesc = row_metadata_merged[, -1, drop = FALSE],
      cid = colnames(norm_counts)
    ),
    ofile = paste0(output_prefix, "_norm_counts.gct")
  )
  
  # Export subset of significant genes if there are any
  if (sum(!is.na(lrt_results$padj) & lrt_results$padj < args$fdr) > 0) {
    log_message("Exporting significant genes count matrix...", "INFO")
    
    # Get significant genes
    sig_genes <- subset_significant_genes(lrt_results, args$fdr)
    
    # Only export if we have significant genes
    if (length(sig_genes) > 0) {
      # Subset normalized counts and metadata
      sig_counts <- norm_counts[sig_genes, , drop = FALSE]
      row_metadata_filtered <- row_metadata_merged[row_metadata_merged$id %in% sig_genes, , drop = FALSE]
      
      # Export GCT
      cmapR::write_gct(
        ds = cmapR::new_gct(
          mat = as.matrix(sig_counts),
          rid = row_metadata_filtered$id,
          rdesc = row_metadata_filtered[, -1, drop = FALSE],
          cid = colnames(sig_counts)
        ),
        ofile = paste0(output_prefix, "_significant_genes.gct")
      )
      
      log_message(glue::glue("Exported GCT with {length(sig_genes)} significant genes"), "SUCCESS")
    } else {
      log_message("No significant genes to export", "WARNING")
    }
  } else {
    log_message("No significant genes found at FDR < {args$fdr}", "WARNING")
  }
}

# Export subset of significant genes only
subset_significant_genes <- function(lrt_results, alpha = 0.1) {
  log_message("=== Debugging before subsetting ===", "DEBUG")
  # Get genes with padj < alpha
  sig_genes <- rownames(lrt_results)[which(lrt_results$padj < alpha)]
  log_message(paste("Checking if all row names in row_metadata_filtered are present in normCounts:"), "DEBUG")
  log_message(paste("All row names present:", !is.null(sig_genes)), "DEBUG")
  log_message(paste("Number of rows after subsetting:", length(sig_genes)), "DEBUG")
  return(sig_genes)
}

# Export analysis parameters
export_parameters <- function(args, rpkm_filtered_count, output_prefix) {
  log_message("Exporting analysis parameters...", "INFO")
  
  # Create parameters list
  params <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    input_files = args$input,
    sample_names = args$name,
    design_formula = args$design,
    reduced_formula = args$reduced,
    rpkm_cutoff = args$rpkm_cutoff,
    min_counts = args$mincounts,
    fdr_threshold = args$fdr,
    batch_correction = args$batchcorrection,
    filtering_stats = list(
      rpkm_filtered_genes = rpkm_filtered_count
    )
  )
  
  # Save as RDS and JSON
  saveRDS(params, file = paste0(output_prefix, "_parameters.rds"))
  
  # Convert to JSON if jsonlite is available
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    writeLines(
      jsonlite::toJSON(params, pretty = TRUE, auto_unbox = TRUE),
      paste0(output_prefix, "_parameters.json")
    )
  }
  
  log_message("Analysis parameters exported successfully", "SUCCESS")
} 