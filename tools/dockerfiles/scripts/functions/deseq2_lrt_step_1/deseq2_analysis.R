#!/usr/bin/env Rscript

# --- DESeq2 analysis functions ---

# Main function to run DESeq2 analysis
run_deseq2 <- function(count_data, metadata, design_formula, reduced_formula, args) {
  logger::info("Running DESeq2 analysis with LRT")
  
  # Create DESeqDataSet object
  dds <- try(DESeqDataSetFromMatrix(
    countData = count_data,
    colData = metadata,
    design = design_formula
  ))
  
  if (inherits(dds, "try-error")) {
    stop("Failed to create DESeqDataSet from matrix")
  }
  
  # Run DESeq with LRT
  logger::info(paste("Running DESeq with full model:", deparse(design_formula)))
  logger::info(paste("Reduced model:", deparse(reduced_formula)))
  
  dds <- try(DESeq(
    dds,
    test = "LRT",
    reduced = reduced_formula,
    parallel = (args$threads > 1),
    BPPARAM = MulticoreParam(args$threads)
  ))
  
  if (inherits(dds, "try-error")) {
    stop("DESeq2 analysis failed")
  }
  
  # Extract results with FDR threshold
  results <- try(results(
    dds,
    alpha = args$fdr,
    lfcThreshold = args$lfcthreshold
  ))
  
  if (inherits(results, "try-error")) {
    stop("Failed to extract results from DESeq2 analysis")
  }
  
  # Sort results by p-value
  results <- results[order(results$pvalue), ]
  
  # Get normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Get VST transformed data
  vst_data <- try(vst(dds, blind = FALSE))
  
  if (inherits(vst_data, "try-error")) {
    logger::warn("VST transformation failed, will attempt to use rlog instead")
    vst_data <- try(rlog(dds, blind = FALSE))
    
    if (inherits(vst_data, "try-error")) {
      logger::error("Both VST and rlog transformations failed")
      vst_data <- NULL
    }
  }
  
  return(list(
    dds = dds,
    results = results,
    normalized_counts = normalized_counts,
    vst_data = vst_data
  ))
}

# Function to export results from DESeq2 analysis
export_results <- function(deseq_results, expression_data, metadata, args) {
  logger::info("Exporting DESeq2 results")
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(args$output_prefix)
  if (!dir.exists(output_dir) && output_dir != "") {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract result components
  dds <- deseq_results$dds
  results <- deseq_results$results
  normalized_counts <- deseq_results$normalized_counts
  vst_data <- deseq_results$vst_data
  
  # Export DESeq2 results table
  results_df <- as.data.frame(results)
  results_df$GeneId <- row.names(results_df)
  
  # Reorder columns to put GeneId first
  results_df <- results_df[, c("GeneId", setdiff(colnames(results_df), "GeneId"))]
  
  # Write results to TSV file
  write.table(
    results_df,
    file = paste0(args$output_prefix, ".deseq2_lrt_results.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Export normalized counts
  norm_counts_df <- as.data.frame(normalized_counts)
  norm_counts_df$GeneId <- row.names(norm_counts_df)
  norm_counts_df <- norm_counts_df[, c("GeneId", setdiff(colnames(norm_counts_df), "GeneId"))]
  
  write.table(
    norm_counts_df,
    file = paste0(args$output_prefix, ".normalized_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Export VST data if available
  if (!is.null(vst_data)) {
    vst_df <- as.data.frame(assay(vst_data))
    vst_df$GeneId <- row.names(vst_df)
    vst_df <- vst_df[, c("GeneId", setdiff(colnames(vst_df), "GeneId"))]
    
    write.table(
      vst_df,
      file = paste0(args$output_prefix, ".vst_transformed.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  # Create basic plots
  # PCA plot
  if (!is.null(vst_data)) {
    pdf(paste0(args$output_prefix, ".pca_plot.pdf"), width = 10, height = 8)
    plotPCA(vst_data, intgroup = c("treatment", "cond"))
    dev.off()
  }
  
  # MA plot
  pdf(paste0(args$output_prefix, ".ma_plot.pdf"), width = 10, height = 8)
  plotMA(results)
  dev.off()
  
  # Volcano plot
  pdf(paste0(args$output_prefix, ".volcano_plot.pdf"), width = 10, height = 8)
  with(results_df, {
    plot(
      log2FoldChange, 
      -log10(pvalue),
      pch = 20,
      main = "Volcano Plot",
      xlab = "log2 Fold Change",
      ylab = "-log10 p-value"
    )
    significant <- padj < args$fdr & abs(log2FoldChange) > args$lfcthreshold
    points(log2FoldChange[significant], -log10(pvalue[significant]), pch = 20, col = "red")
  })
  dev.off()
  
  logger::info("Results exported successfully")
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

# Function to validate essential arguments
validate_essential_args <- function(args) {
  # Check if args is a valid list
  if (!is.list(args)) {
    stop("Arguments must be provided as a list")
  }
  
  # Check essential arguments
  essential_args <- c("meta", "design", "reduced", "input", "name")
  missing_args <- essential_args[!essential_args %in% names(args)]
  
  if (length(missing_args) > 0) {
    stop(paste("Missing essential arguments:", paste(missing_args, collapse=", ")))
  }
  
  # Check if any essential arguments are NULL or empty
  empty_args <- c()
  for (arg in essential_args) {
    if (is.null(args[[arg]]) || (is.character(args[[arg]]) && length(args[[arg]]) == 0)) {
      empty_args <- c(empty_args, arg)
    }
  }
  
  if (length(empty_args) > 0) {
    stop(paste("These essential arguments are empty:", paste(empty_args, collapse=", ")))
  }
  
  # Check that input and name have equal lengths
  if (length(args$input) != length(args$name)) {
    logger::warn(paste("Mismatch between number of input files (", length(args$input), 
                        ") and sample names (", length(args$name), ")"))
    
    # If more samples than inputs, trim sample names
    if (length(args$name) > length(args$input)) {
      logger::warn("Trimming extra sample names to match input files")
      args$name <- args$name[1:length(args$input)]
    }
    
    # If more inputs than samples, trim input files
    if (length(args$input) > length(args$name)) {
      logger::warn("Trimming extra input files to match sample names")
      args$input <- args$input[1:length(args$name)]
    }
  }
  
  logger::info("Essential arguments validated successfully")
  return(args)
}

# Main function to perform DESeq2 LRT analysis
run_deseq2_lrt_analysis <- function(args) {
  # Log the start of the analysis
  logger::info("Starting DESeq2 analysis with LRT")
  logger::debug("Arguments:", args)
  
  # Validate essential arguments
  validate_essential_args(args)
  
  # Load the metadata file
  metadata <- load_metadata(args$meta)
  
  # Load expression data
  if (is.null(args$input) || length(args$input) == 0) {
    stop("No input expression files provided")
  }
  
  if (is.null(args$name) || length(args$name) == 0) {
    stop("No sample names provided")
  }
  
  # Make sure name doesn't contain non-sample names
  if (length(args$name) > length(args$input)) {
    logger::warn("More sample names than input files provided. Trimming extra names.")
    args$name <- args$name[1:length(args$input)]
  }
  
  # Now proceed with loading the data
  expression_data <- load_expression_data(args$input, args$name)
  
  # Extract read counts
  counts_cols <- grep(READ_COL, colnames(expression_data), ignore.case = TRUE, value = TRUE)
  if (length(counts_cols) == 0) {
    stop("No count columns found in expression data")
  }
  
  count_data <- expression_data[, c(INTERSECT_BY, counts_cols)]
  
  # Set gene IDs as row names
  row.names(count_data) <- count_data[[INTERSECT_BY]]
  count_data[[INTERSECT_BY]] <- NULL
  
  # Check if there are any zero rows
  zero_rows <- rowSums(count_data) == 0
  if (any(zero_rows)) {
    logger::info(paste("Removing", sum(zero_rows), "rows with zero counts"))
    count_data <- count_data[!zero_rows, , drop = FALSE]
  }
  
  # Set up DESeq2 design from formula string
  if (is.null(args$design) || args$design == "") {
    stop("No valid design formula provided")
  }
  design_formula <- as.formula(args$design)
  
  # Set up reduced design for LRT
  if (is.null(args$reduced) || args$reduced == "") {
    stop("No valid reduced design formula provided")
  }
  reduced_formula <- as.formula(args$reduced)
  
  # Prepare metadata and ensure it matches count data samples
  metadata <- prepare_metadata_for_deseq(metadata, colnames(count_data))
  
  # Run DESeq2 analysis
  deseq_results <- run_deseq2(count_data, metadata, design_formula, reduced_formula, args)
  
  # Export results
  export_results(deseq_results, expression_data, metadata, args)
} 