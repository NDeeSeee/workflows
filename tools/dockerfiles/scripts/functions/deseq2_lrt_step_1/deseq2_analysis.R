#!/usr/bin/env Rscript

# --- DESeq2 analysis functions ---

# Main function to run DESeq2 analysis
run_deseq2 <- function(count_data, sample_metadata, design_formula, reduced_formula, args) {
  message("Running DESeq2 analysis...")
  
  # Set up parallel processing if requested
  if (args$threads > 1) {
    message(paste("Using", args$threads, "threads for parallel processing"))
    BiocParallel::register(BiocParallel::MulticoreParam(args$threads))
  }
  
  # Create DESeqDataSet
  message("Creating DESeqDataSet object")
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_metadata,
    design = design_formula
  )
  
  # Filter low count genes if requested
  if (!is.null(args$filter_low_counts) && args$filter_low_counts > 0) {
    message(paste("Filtering genes with fewer than", args$filter_low_counts, "counts"))
    keep <- rowSums(counts(dds)) >= args$filter_low_counts
    dds <- dds[keep,]
    message(paste("Retained", sum(keep), "of", length(keep), "genes after filtering"))
  }
  
  # Run DESeq2 with LRT test
  message("Running DESeq with LRT test")
  message(paste("Design formula:", deparse(design_formula)))
  message(paste("Reduced formula:", deparse(as.formula(reduced_formula))))
  
  dds <- DESeq(
    dds,
    test = "LRT",
    reduced = as.formula(reduced_formula),
    parallel = (args$threads > 1),
    quiet = FALSE
  )
  
  # Extract results
  message("Extracting results")
  res <- results(
    dds,
    alpha = args$fdr,
    pAdjustMethod = "BH"
  )
  
  # Sort results by p-value
  ordered_res <- res[order(res$pvalue),]
  
  # Annotate with normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Apply transformation based on data size
  if (ncol(dds) >= 30) {
    # Use VST for larger datasets (faster)
    message("Applying variance stabilizing transformation (VST)")
    vst_data <- vst(dds, blind = FALSE)
    transformed_counts <- assay(vst_data)
  } else {
    # Use rlog for smaller datasets (more accurate)
    message("Applying regularized log transformation (rlog)")
    rlog_data <- rlog(dds, blind = FALSE)
    transformed_counts <- assay(rlog_data)
  }
  
  # Return results
  return(list(
    dds = dds,
    res = res,
    ordered_res = ordered_res,
    normalized_counts = normalized_counts,
    transformed_counts = transformed_counts
  ))
}

# Function to export DESeq2 results
export_results <- function(deseq_results, expression_df, metadata_df, args, batch_warning = NULL, filtered_count = NULL) {
  message("Exporting DESeq2 results...")
  
  # Create output directory if it doesn't exist
  output_dir <- args$output_prefix
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Get results components
  dds <- deseq_results$dds
  res <- deseq_results$res
  ordered_res <- deseq_results$ordered_res
  normalized_counts <- deseq_results$normalized_counts
  transformed_counts <- deseq_results$transformed_counts
  
  # 1. Export detailed results table with annotations
  message("Exporting detailed DESeq2 results table")
  
  # Convert results to data frame and add gene information
  results_df <- as.data.frame(ordered_res)
  results_df$GeneId <- rownames(results_df)
  
  # Add gene annotations from expression data if available
  if ("GeneSymbol" %in% colnames(expression_df)) {
    gene_info <- expression_df[, c("GeneId", "GeneSymbol"), drop = FALSE]
    gene_info <- gene_info[!duplicated(gene_info$GeneId), ]
    results_df <- merge(results_df, gene_info, by = "GeneId", all.x = TRUE, sort = FALSE)
  }
  
  # Add log10 transformed p-values for better visualization
  results_df$log10_pvalue <- -log10(results_df$pvalue)
  results_df$log10_padj <- -log10(results_df$padj)
  
  # Reorder columns for better readability
  if ("GeneSymbol" %in% colnames(results_df)) {
    col_order <- c("GeneId", "GeneSymbol", "baseMean", "log2FoldChange", "lfcSE", 
                  "stat", "pvalue", "padj", "log10_pvalue", "log10_padj")
  } else {
    col_order <- c("GeneId", "baseMean", "log2FoldChange", "lfcSE", 
                  "stat", "pvalue", "padj", "log10_pvalue", "log10_padj")
  }
  
  # Reorder columns, keeping any additional columns at the end
  other_cols <- base::setdiff(colnames(results_df), col_order)
  col_order <- c(col_order, other_cols)
  results_df <- results_df[, base::intersect(col_order, colnames(results_df)), drop = FALSE]
  
  # Write results to TSV file
  results_file <- file.path(output_dir, "deseq2_lrt_results.tsv")
  write.table(results_df, file = results_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message(paste("Wrote results to:", results_file))
  
  # 2. Export normalized counts
  message("Exporting normalized counts")
  norm_counts_df <- as.data.frame(normalized_counts)
  norm_counts_df$GeneId <- rownames(norm_counts_df)
  
  # Move GeneId to first column
  norm_counts_df <- norm_counts_df[, c("GeneId", base::setdiff(colnames(norm_counts_df), "GeneId")), drop = FALSE]
  
  # Write normalized counts to TSV file
  norm_counts_file <- file.path(output_dir, "normalized_counts.tsv")
  write.table(norm_counts_df, file = norm_counts_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message(paste("Wrote normalized counts to:", norm_counts_file))
  
  # 3. Export transformed counts
  message("Exporting transformed counts")
  trans_counts_df <- as.data.frame(transformed_counts)
  trans_counts_df$GeneId <- rownames(trans_counts_df)
  
  # Move GeneId to first column
  trans_counts_df <- trans_counts_df[, c("GeneId", base::setdiff(colnames(trans_counts_df), "GeneId")), drop = FALSE]
  
  # Write transformed counts to TSV file
  trans_counts_file <- file.path(output_dir, "transformed_counts.tsv")
  write.table(trans_counts_df, file = trans_counts_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message(paste("Wrote transformed counts to:", trans_counts_file))
  
  # 4. Export GCT file for downstream analysis
  message("Exporting GCT file for downstream analysis")
  
  # Create matrix for GCT export
  gct_data <- transformed_counts
  
  # Safely extract arguments for scaling and clustering
  scaling_type <- if (!is.null(args$scaling_type)) as.character(args$scaling_type) else "none"
  cluster_method <- if (!is.null(args$cluster_method)) as.character(args$cluster_method) else "none"
  row_distance <- if (!is.null(args$row_distance)) as.character(args$row_distance) else "euclid"
  column_distance <- if (!is.null(args$column_distance)) as.character(args$column_distance) else "euclid"
  k_hopach <- if (!is.null(args$k_hopach)) as.numeric(args$k_hopach) else 3
  kmax_hopach <- if (!is.null(args$kmax_hopach)) as.numeric(args$kmax_hopach) else 5
  
  # Validate parameters
  if (!scaling_type %in% c("none", "zscore", "minmax")) {
    warning(paste("Unknown scaling type:", scaling_type, "- defaulting to none"))
    scaling_type <- "none"
  }
  
  if (!cluster_method %in% c("none", "row", "column", "both")) {
    warning(paste("Unknown cluster method:", cluster_method, "- defaulting to none"))
    cluster_method <- "none"
  }
  
  # Scale data if requested
  if (scaling_type != "none") {
    tryCatch({
      gct_data <- scale_expression_data(gct_data, scaling_type)
    }, error = function(e) {
      warning(paste("Error during scaling:", e$message, "- using unscaled data"))
    })
  }
  
  # Apply clustering if requested
  if (cluster_method != "none") {
    tryCatch({
      clustering_result <- perform_clustering(
        gct_data,
        cluster_method,
        row_distance,
        column_distance,
        k_hopach,
        kmax_hopach
      )
      gct_data <- clustering_result$expression_data
    }, error = function(e) {
      warning(paste("Error during clustering:", e$message, "- using unclustered data"))
    })
  }
  
  # Write GCT file
  gct_file <- file.path(output_dir, "deseq2_transformed_data.gct")
  write_gct(gct_data, gct_file)
  message(paste("Wrote GCT file to:", gct_file))
  
  # 5. Export summary report
  message("Generating summary report")
  summary_file <- file.path(output_dir, "analysis_summary.txt")
  
  # Count significant genes
  sig_genes <- sum(!is.na(res$padj) & res$padj < args$fdr, na.rm = TRUE)
  
  # Write summary to file
  writeLines(
    c(
      "DESeq2 LRT Analysis Summary",
      "==========================",
      "",
      paste("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("Design Formula:", deparse(design(dds))),
      paste("Reduced Formula:", args$reduced),
      paste("FDR Threshold:", args$fdr),
      paste("LFC Threshold:", args$lfcthreshold),
      "",
      paste("Total Genes Analyzed:", nrow(dds)),
      paste("Significant Genes (FDR <", args$fdr, "):", sig_genes),
      paste("Filtering:", if(!is.null(filtered_count)) paste(filtered_count, "genes removed by RPKM filtering") else "No RPKM filtering applied"),
      paste("Batch Correction:", if(!is.null(args$batchcorrection)) args$batchcorrection else "none"),
      if(!is.null(batch_warning)) c("", paste("Warning:", batch_warning)) else character(0)
    ),
    con = summary_file
  )
  message(paste("Wrote summary report to:", summary_file))
  
  message("Results export completed")
}

# Function to run DESeq2 LRT analysis
run_deseq2_lrt <- function(count_data, sample_metadata, design_formula, reduced_formula, 
                           batchcorrection = "none", threads = 1) {
  require(DESeq2)
  
  # Create DESeqDataSet
  message("Creating DESeqDataSet...")
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_metadata,
    design = design_formula
  )
  
  # Apply batch correction if requested
  if (batchcorrection == "combatseq") {
    message("Applying ComBat-Seq batch correction...")
    require(sva)
    
    # Check if 'batch' column exists in sample metadata
    if (!("batch" %in% colnames(sample_metadata))) {
      stop("Batch correction requested but 'batch' column not found in metadata")
    }
    
    # Get batch information
    batch <- sample_metadata$batch
    
    # Get model matrix from full design
    group <- sample_metadata[, sapply(sample_metadata, is.factor)]
    mod <- model.matrix(design_formula, data = group)
    
    # Apply ComBat-Seq
    corrected_counts <- ComBat_seq(
      counts = count_data,
      batch = batch,
      group = group
    )
    
    # Update DESeqDataSet with corrected counts
    dds <- DESeqDataSetFromMatrix(
      countData = corrected_counts,
      colData = sample_metadata,
      design = design_formula
    )
  }
  
  # Set up parallel processing
  if (threads > 1) {
    message(paste("Using", threads, "threads for parallel processing"))
    BiocParallel::register(BiocParallel::MulticoreParam(threads))
  }
  
  # Run DESeq2 with LRT
  message("Running DESeq2 LRT analysis...")
  dds <- DESeq(
    dds,
    test = "LRT",
    reduced = reduced_formula,
    parallel = (threads > 1)
  )
  
  return(dds)
}

# Function to extract results from DESeq2 analysis
extract_deseq2_results <- function(dds, fdr = 0.1, lfcthreshold = 0.59, use_lfc_thresh = FALSE) {
  require(DESeq2)
  
  # Determine if we're dealing with LRT or Wald test
  is_lrt <- mcols(dds)$modelMatrixType == "LRT"
  
  if (is_lrt) {
    message(paste("Extracting DESeq2 LRT results with FDR =", fdr))
    
    # For LRT test, don't use lfcThreshold parameter
    res <- results(
      dds,
      alpha = fdr,
      altHypothesis = "greaterAbs",
      pAdjustMethod = "BH",
      independentFiltering = TRUE,
      cooksCutoff = TRUE,
      parallel = FALSE
    )
  } else {
    # Wald test logic
    if (use_lfc_thresh) {
      message(paste("Extracting DESeq2 Wald results with FDR =", fdr, 
                    "and LFC threshold =", lfcthreshold, "as hypothesis threshold"))
      
      # If use_lfc_thresh is TRUE, use lfcthreshold directly in hypothesis test
      res <- results(
        dds,
        alpha = fdr,
        lfcThreshold = lfcthreshold,
        altHypothesis = "greaterAbs",
        pAdjustMethod = "BH",
        independentFiltering = TRUE,
        cooksCutoff = TRUE,
        parallel = FALSE
      )
    } else {
      message(paste("Extracting DESeq2 Wald results with FDR =", fdr, 
                    "and post-filtering with LFC threshold =", lfcthreshold))
      
      # If use_lfc_thresh is FALSE, use lfcThreshold = 0 and filter later
      res <- results(
        dds,
        alpha = fdr,
        lfcThreshold = 0,
        altHypothesis = "greaterAbs",
        pAdjustMethod = "BH",
        independentFiltering = TRUE,
        cooksCutoff = TRUE,
        parallel = FALSE
      )
    }
  }
  
  # Convert to data frame for easier manipulation
  res_df <- as.data.frame(res)
  
  # Add gene names if available
  if (!is.null(mcols(dds)$symbol)) {
    res_df$symbol <- mcols(dds)$symbol
  }
  
  # Apply post-filtering for LFC if use_lfc_thresh is FALSE and it's not an LRT test
  if (!is_lrt && !use_lfc_thresh && lfcthreshold > 0) {
    message(paste("Post-filtering results with |log2FoldChange| >=", lfcthreshold))
    
    # Mark genes not meeting threshold as non-significant
    significant_genes <- !is.na(res_df$padj) & res_df$padj < fdr & abs(res_df$log2FoldChange) >= lfcthreshold
    
    # Create a new column for LFC-filtered p-values
    res_df$padj_lfc_filtered <- res_df$padj
    res_df$padj_lfc_filtered[!significant_genes] <- NA
    
    # Create a filtering status column
    res_df$lfc_significant <- significant_genes
    
    message(paste("After LFC filtering:", sum(significant_genes, na.rm = TRUE), "of", 
                  sum(!is.na(res_df$padj) & res_df$padj < fdr), "significant genes remain"))
  }
  
  # Sort by adjusted p-value
  res_df <- res_df[order(res_df$padj), ]
  
  return(res_df)
}

# Function to apply limma batch correction on the results
apply_limma_batch_correction <- function(dds, sample_metadata) {
  require(limma)
  
  message("Applying limma batch correction to normalized counts...")
  
  # Check if 'batch' column exists in sample metadata
  if (!("batch" %in% colnames(sample_metadata))) {
    stop("Batch correction requested but 'batch' column not found in metadata")
  }
  
  # Get normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Create design matrix
  model <- model.matrix(~ 0 + group, data = sample_metadata)
  colnames(model) <- levels(sample_metadata$group)
  
  # Create batch variable
  batch <- sample_metadata$batch
  
  # Apply removeBatchEffect
  corrected_counts <- removeBatchEffect(
    normalized_counts,
    batch = batch,
    design = model
  )
  
  return(corrected_counts)
} 