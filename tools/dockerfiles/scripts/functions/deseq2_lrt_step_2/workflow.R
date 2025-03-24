#!/usr/bin/env Rscript
#
# Workflow functions for DESeq2 LRT Step 2
#

#' Main workflow function for DESeq2 LRT Step 2
#'
#' @param args Command-line arguments
#' @return Final results list containing normalized counts, expression data, and contrast results
#' @export
run_workflow <- function(args) {
  log_message("Starting DESeq2 LRT Step 2 workflow")
  
  # Load input data from Step 1
  log_message(paste("Loading DESeq2 object data from", args$dsq_obj_data))
  step1_data <- readRDS(args$dsq_obj_data)
  
  # Load contrast data
  log_message(paste("Loading contrast data from", args$contrast_df))
  contrast_df <- read.delim(args$contrast_df, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Get contrast indices to process
  contrast_indices <- as.numeric(unlist(strsplit(args$contrast_indices, ",")))
  log_message(paste("Processing contrast indices:", paste(contrast_indices, collapse = ", ")))
  
  # Check if contrast indices are valid
  if (any(contrast_indices > nrow(contrast_df))) {
    stop("Invalid contrast index. Maximum index is ", nrow(contrast_df))
  }
  
  # Subset contrast data to process only the specified contrasts
  contrast_df <- contrast_df[contrast_indices, , drop = FALSE]
  
  # Initialize variables
  dds <- step1_data$dds
  expDataDf <- step1_data$expDataDf
  
  # Check for required columns in the DESeq dataset
  if (is.null(dds)) {
    stop("DESeq2 dataset is missing from step 1 data")
  }
  
  # Extract metadata from DESeq2 object
  col_data <- as.data.frame(colData(dds))
  design_formula <- design(dds)
  
  # Verify metadata and consistency between steps
  col_data <- validate_metadata(col_data, args$batchcorrection, design_formula)
  verify_step_consistency(col_data, design_formula, args)
  
  # Identify all factors in the DESeq dataset for reference
  factor_cols <- sapply(col_data, is.factor)
  factor_names <- names(factor_cols)[factor_cols]
  log_message(paste("Available factors in DESeq2 dataset:", paste(factor_names, collapse = ", ")))
  
  # Rebuild dds with reference levels if needed
  reference_levels <- list()
  if (exists("factors", where = step1_data)) {
    for (factor_name in step1_data$factors) {
      if (factor_name %in% factor_names) {
        reference_levels[[factor_name]] <- levels(dds[[factor_name]])[1]
      } else {
        log_warning(paste("Factor", factor_name, "from step 1 not found in DESeq2 dataset"))
      }
    }
    log_message(paste("Using reference levels:", paste(names(reference_levels), reference_levels, sep = "=", collapse = ", ")))
  } else {
    log_warning("No factors list found in step 1 data - using factor information from DESeq2 dataset")
    for (factor_name in factor_names) {
      if (factor_name != "batch") {  # Skip batch factor
        reference_levels[[factor_name]] <- levels(dds[[factor_name]])[1]
      }
    }
  }
  
  # Get normalized counts from DESeq2
  normCounts <- counts(dds, normalized = TRUE)
  log_message(paste("Extracted normalized counts with dimensions:", nrow(normCounts), "x", ncol(normCounts)))
  
  # Apply batch correction if requested
  if (args$batchcorrection != "none") {
    log_message(paste("Applying", args$batchcorrection, "batch correction to normalized counts"))
    normCounts <- apply_batch_correction(
      count_data = normCounts,
      metadata_df = col_data,
      batch_method = args$batchcorrection,
      normalized = TRUE  # These are normalized counts
    )
    log_message("Batch correction completed")
  } else {
    log_message("No batch correction applied to normalized counts")
  }
  
  # Create MDS plot for sample visualization
  log_message("Creating MDS plot for sample visualization")
  with_error_handling({
    create_mds_plot(
      normCounts = normCounts,
      col_metadata = col_data,
      output_file = "mds_plot.html",
      interactive = TRUE
    )
  })
  
  # Process each contrast
  results_list <- list()
  
  for (i in 1:nrow(contrast_df)) {
    contrast_row <- contrast_df[i, ]
    contrast_name <- contrast_row$contrast_name
    log_message(paste("Processing contrast", i, "of", nrow(contrast_df), ":", contrast_name))
    
    # Get contrast results
    contrast_res <- with_error_handling({
      get_contrast_res(step1_data, contrast_row, args)
    })
    
    if (is.null(contrast_res)) {
      log_warning(paste("Failed to process contrast:", contrast_name))
      next
    }
    
    # Add results to expression data
    updated_expDataDf <- with_error_handling({
      add_metadata_to_results(expDataDf, contrast_res, contrast_name)
    })
    
    if (!is.null(updated_expDataDf)) {
      expDataDf <- updated_expDataDf
      log_message(paste("Added results for contrast", contrast_name, "to expression data"))
    }
    
    # Store results
    results_list[[contrast_name]] <- contrast_res
    
    # Force garbage collection after each contrast to manage memory
    gc(verbose = FALSE)
  }
  
  # Prepare final results
  final_results <- list(
    dds = dds,
    normCounts = normCounts,
    expDataDf = expDataDf,
    results = results_list,
    contrasts = contrast_df
  )
  
  # Save results
  if (!is.null(args$output)) {
    output_file <- paste0(args$output, "_deseq2_step2_results.rds")
    log_message(paste("Saving results to", output_file))
    saveRDS(final_results, file = output_file)
  }
  
  # Export reports
  log_message("Exporting reports")
  
  # Export gene expression table
  export_deseq_report(expDataDf, args$output)
  
  # Export GCT data for visualization
  row_metadata <- expDataDf
  # Ensure gene_id is used as row names if present
  if ("gene_id" %in% colnames(row_metadata)) {
    rownames(row_metadata) <- row_metadata$gene_id
  } else if ("GeneId" %in% colnames(row_metadata)) {
    rownames(row_metadata) <- row_metadata$GeneId
  }
  
  col_metadata <- col_data
  
  export_gct_data(normCounts, row_metadata, col_metadata, args)
  
  log_message("DESeq2 LRT Step 2 workflow completed successfully")
  
  return(final_results)
} 