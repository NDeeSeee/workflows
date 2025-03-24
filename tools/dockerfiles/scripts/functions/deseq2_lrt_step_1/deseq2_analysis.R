#!/usr/bin/env Rscript

# --- DESeq2 analysis functions ---

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