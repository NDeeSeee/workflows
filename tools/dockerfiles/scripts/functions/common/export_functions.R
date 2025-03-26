#!/usr/bin/env Rscript
#
# Consolidated export functions for DESeq2 LRT analysis
#

#' Export DESeq2 report as a TSV file
#'
#' @param expression_data_df Expression data with DESeq2 results
#' @param output_prefix Prefix for output filename
#' @export
export_deseq_report <- function(expression_data_df, output_prefix) {
  expression_data_df_filename <- paste0(output_prefix, "_gene_exp_table.tsv")
  write.table(
    expression_data_df,
    file = expression_data_df_filename,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  log_message(paste("Exported DESeq report to", expression_data_df_filename))
}

#' Export CLS file for GSEA
#'
#' @param categories Factor with category labels
#' @param location Output file location
#' @export
export_cls <- function(categories, location) {
  tryCatch(
    expr = {
      output_stream <- file(location, "w")
      on.exit(close(output_stream), add = TRUE)
      cat(
        paste(length(categories),
              length(levels(categories)),
              "1",
              sep = "\t"
        ),
        paste("#", paste(
          unique(as.character(categories)),
          collapse = "\t"
        ), sep = "\t"),
        paste(paste(as.character(categories),
                    collapse = "\t"
        ), sep = "\t"),
        file = output_stream,
        sep = "\n"
      )
      log_message(paste("Exported CLS data to", location))
    },
    error = function(e) {
      log_error(paste("Failed to export CLS data to", location, "with error -", e$message))
    }
  )
}

#' Export results and visualizations for Step 1
#'
#' @param deseq_results Results from DESeq2 analysis
#' @param expression_df Expression data frame
#' @param metadata_df Metadata data frame
#' @param args Command-line arguments
#' @param batch_warning Any warnings from batch correction
#' @param rpkm_filtered_count Count of genes filtered by RPKM
#' @export
export_results <- function(deseq_results, expression_df, metadata_df, args, batch_warning, rpkm_filtered_count) {
  log_message("Exporting results and visualizations...", "STEP")
  
  # No need to harmonize parameters as they're already handled in the ArgumentParser
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(args$output_prefix)
  if (!dir.exists(output_dir) && output_dir != ".") {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set output prefix
  output_prefix <- args$output_prefix
  
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

#' Export data table as TSV
#'
#' @param lrt_res LRT results from DESeq2
#' @param expression_df Expression data frame
#' @param output_prefix Output file prefix
#' @param args Command-line arguments
#' @export
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

#' Export contrast information
#'
#' @param contrasts Contrast data frame
#' @param output_prefix Output file prefix
#' @export
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

#' Export visualizations
#'
#' @param dds DESeq2 dataset
#' @param lrt_results LRT results from DESeq2
#' @param metadata_df Metadata data frame
#' @param norm_counts Normalized counts matrix
#' @param output_prefix Output file prefix
#' @param args Command-line arguments
#' @export
export_visualizations <- function(dds, lrt_results, metadata_df, norm_counts, output_prefix, args) {
  log_message("Generating visualizations...", "STEP")
  
  # Generate MDS plot
  mds_plot_file <- paste0(output_prefix, "_mds_plot.html")
  create_mds_plot(
    normCounts = norm_counts,
    col_metadata = metadata_df,
    output_file = mds_plot_file,
    interactive = TRUE
  )
  
  # Generate MA plot
  if (!is.null(lrt_results)) {
    ma_plot <- generate_ma_plot(lrt_results)
    ggsave(paste0(output_prefix, "_ma_plot.png"), plot = ma_plot, width = 10, height = 8)
    
    # Generate p-value distribution
    pval_plot <- generate_pvalue_histogram(lrt_results)
    ggsave(paste0(output_prefix, "_pvalue_hist.png"), plot = pval_plot, width = 10, height = 6)
  }
  
  log_message("Visualizations exported successfully", "SUCCESS")
}

#' Export normalized counts and filtered counts to GCT format
#'
#' @param normCounts Normalized counts matrix
#' @param row_metadata Row metadata data frame
#' @param col_metadata Column metadata data frame
#' @param args Command-line arguments or output prefix
#' @param output_prefix Optional output prefix (used if args doesn't contain output prefix)
#' @export
export_gct_data <- function(normCounts, row_metadata, col_metadata, args, output_prefix = NULL) {
  tryCatch({
    log_message("Exporting GCT data...", "STEP")
    
    # Determine filename based on args or output_prefix
    if (is.character(args)) {
      output_prefix <- args
      args <- list(fdr = 0.1, lfcthreshold = 1)
    } else if (is.null(output_prefix)) {
      # Check for both old and new parameter names
      if (!is.null(args$output_prefix)) {
        output_prefix <- args$output_prefix
      } else if (!is.null(args$output)) {
        output_prefix <- args$output
      }
    }
    
    # Default filenames
    all_counts_file <- "counts_all.gct"
    filtered_counts_file <- "counts_filtered.gct"
    
    # Use prefix if provided
    if (!is.null(output_prefix)) {
      all_counts_file <- paste0(output_prefix, "_", all_counts_file)
      filtered_counts_file <- paste0(output_prefix, "_", filtered_counts_file)
    }
    
    # Ensure col_metadata columns are vectors
    if (!is.null(col_metadata) && nrow(col_metadata) > 0) {
      col_metadata <- col_metadata %>% dplyr::mutate_all(as.vector)
      
      # Ensure sample ordering is consistent
      common_cols <- intersect(colnames(normCounts), rownames(col_metadata))
      if (length(common_cols) > 0) {
        normCounts <- normCounts[, common_cols, drop = FALSE]
        col_metadata <- col_metadata[common_cols, , drop = FALSE]
      }
    }
    
    # Check dimensions
    if (nrow(normCounts) == 0 || ncol(normCounts) == 0) {
      log_error("Cannot export empty count matrix")
      return(NULL)
    }
    
    # Fix rownames in metadata if needed
    if ("gene_id" %in% colnames(row_metadata)) {
      rownames(row_metadata) <- row_metadata$gene_id
    } else if ("GeneId" %in% colnames(row_metadata)) {
      rownames(row_metadata) <- row_metadata$GeneId
    }
    
    # Create common row indices between counts and metadata
    common_rows <- intersect(rownames(normCounts), rownames(row_metadata))
    if (length(common_rows) == 0) {
      log_error("No common rows between counts and metadata")
      return(NULL)
    }
    
    # Subset data to common rows
    normCounts_common <- normCounts[common_rows, , drop = FALSE]
    row_metadata_common <- row_metadata[common_rows, , drop = FALSE]
    
    # Export all counts
    gct_data <- cmapR::GCT(
      mat = as.matrix(normCounts_common),
      rdesc = as.data.frame(row_metadata_common),
      cdesc = as.data.frame(col_metadata)
    )
    cmapR::write_gct(gct_data, all_counts_file)
    log_message(paste("Exported GCT (all) to", all_counts_file))
    
    # Filter rows by any FDR column
    fdr_cols <- grep("_FDR$|padj", colnames(row_metadata), value = TRUE)
    lfc_cols <- grep("_LFC$|log2FoldChange", colnames(row_metadata), value = TRUE)
    
    # Apply filtering if FDR columns exist
    if (length(fdr_cols) == 0) {
      log_warning("No FDR columns found. No filtering performed for filtered GCT.")
      return(all_counts_file)
    } 
    
    # Filter by FDR and LFC
    row_metadata_filtered <- row_metadata %>%
      dplyr::filter(dplyr::if_any(dplyr::all_of(fdr_cols), ~ !is.na(.) & . <= args$fdr))
    
    # Apply LFC filter if LFC columns exist
    if (length(lfc_cols) > 0) {
      row_metadata_filtered <- row_metadata_filtered %>%
        dplyr::filter(dplyr::if_any(dplyr::all_of(lfc_cols), ~ !is.na(.) & abs(.) >= args$lfcthreshold))
    }
    
    # Check if any genes passed filtering
    if (nrow(row_metadata_filtered) == 0) {
      log_warning(paste("No genes passed the FDR/LFC thresholds of", args$fdr, "/", args$lfcthreshold))
      return(all_counts_file)
    }
    
    log_message(paste(nrow(row_metadata_filtered), "genes passed filtering criteria"))
    
    # Get filtered counts
    filtered_rows <- intersect(rownames(normCounts), rownames(row_metadata_filtered))
    filtered_normCounts <- normCounts[filtered_rows, , drop = FALSE]
    filtered_row_metadata <- row_metadata_filtered[filtered_rows, , drop = FALSE]
    
    # Cluster filtered data if possible
    if (exists("cluster_and_reorder") && nrow(filtered_normCounts) > 1 && ncol(filtered_normCounts) > 1) {
      # Ensure both old and new parameter names are available
      if (!is.null(args$cluster_method)) {
        cluster_param <- args$cluster_method
      } else if (!is.null(args$cluster)) {
        cluster_param <- args$cluster
      } else {
        cluster_param <- "none"
      }
      
      # Only cluster if requested
      if (cluster_param != "none") {
        clustered_data <- cluster_and_reorder(filtered_normCounts, col_metadata, filtered_row_metadata, args)
        
        gct_data_filtered <- cmapR::GCT(
          mat = as.matrix(clustered_data$normCounts),
          rdesc = as.data.frame(clustered_data$row_metadata),
          cdesc = as.data.frame(clustered_data$col_metadata)
        )
      } else {
        gct_data_filtered <- cmapR::GCT(
          mat = as.matrix(filtered_normCounts),
          rdesc = as.data.frame(filtered_row_metadata),
          cdesc = as.data.frame(col_metadata)
        )
      }
    } else {
      gct_data_filtered <- cmapR::GCT(
        mat = as.matrix(filtered_normCounts),
        rdesc = as.data.frame(filtered_row_metadata),
        cdesc = as.data.frame(col_metadata)
      )
    }
    
    cmapR::write_gct(gct_data_filtered, filtered_counts_file)
    log_message(paste("Exported GCT (filtered) to", filtered_counts_file))
    
    return(c(all_counts_file, filtered_counts_file))
  }, error = function(e) {
    log_error(paste("Failed to export GCT data:", e$message))
    return(NULL)
  })
}

#' Export analysis parameters
#'
#' @param args Command-line arguments
#' @param rpkm_filtered_count Count of genes filtered by RPKM
#' @param output_prefix Output file prefix
#' @export
export_parameters <- function(args, rpkm_filtered_count, output_prefix) {
  log_message("Exporting analysis parameters...", "INFO")
  
  # Create parameters list
  params <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    input_files = if ("input" %in% names(args)) args$input else NULL,
    sample_names = if ("name" %in% names(args)) args$name else NULL,
    design_formula = if ("design" %in% names(args)) args$design else NULL,
    reduced_formula = if ("reduced" %in% names(args)) args$reduced else NULL,
    rpkm_cutoff = if ("rpkm_cutoff" %in% names(args)) args$rpkm_cutoff else NULL,
    min_counts = if ("mincounts" %in% names(args)) args$mincounts else NULL,
    fdr_threshold = if ("fdr" %in% names(args)) args$fdr else 0.1,
    lfc_threshold = if ("lfcthreshold" %in% names(args)) args$lfcthreshold else 1,
    batch_correction = if ("batchcorrection" %in% names(args)) args$batchcorrection else NULL,
    clustering_method = if ("cluster_method" %in% names(args)) args$cluster_method else 
                        if ("cluster" %in% names(args)) args$cluster else "none",
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

#' Generate MA plot for DESeq2 results
#'
#' @param lrt_results LRT results from DESeq2
#' @return ggplot object
#' @export
generate_ma_plot <- function(lrt_results) {
  df <- as.data.frame(lrt_results)
  df$significant <- ifelse(!is.na(df$padj) & df$padj < 0.1, "Significant", "Not significant")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log10(baseMean), y = log2FoldChange, color = significant)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_manual(values = c("Significant" = "red", "Not significant" = "black")) +
    ggplot2::labs(
      title = "MA Plot",
      x = "log10(Mean of normalized counts)",
      y = "log2 Fold Change"
    ) +
    ggplot2::theme_bw()
  
  return(p)
}

#' Generate p-value histogram for DESeq2 results
#'
#' @param lrt_results LRT results from DESeq2
#' @return ggplot object
#' @export
generate_pvalue_histogram <- function(lrt_results) {
  df <- as.data.frame(lrt_results)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    ggplot2::labs(
      title = "P-value Distribution",
      x = "P-value",
      y = "Count"
    ) +
    ggplot2::theme_bw()
  
  return(p)
} 