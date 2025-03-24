#!/usr/bin/env Rscript

# --- Main workflow functions ---

# Load and validate metadata
load_and_validate_metadata <- function(args) {
  message("Loading metadata...")
  
  # Check if metadata file is formatted correctly
  check_file_delimiter(args$meta)
  
  # Load metadata
  metadata_df <- read.table(
    args$meta,
    sep = get_file_type(args$meta),
    header = TRUE,
    stringsAsFactors = FALSE,
    row.names = 1
  )
  
  # Clean metadata column and row names
  colnames(metadata_df) <- clean_sample_names(colnames(metadata_df))
  rownames(metadata_df) <- clean_sample_names(rownames(metadata_df))
  
  message(glue::glue("Loaded metadata for {nrow(metadata_df)} samples with {ncol(metadata_df)} covariates"))
  
  # Validate metadata
  check_batch_column(metadata_df, args$batchcorrection)
  
  # Check design formulas
  design_formula <- check_design_formula(args$design, "Full Design-Formula")
  reduced_formula <- check_design_formula(args$reduced, "Reduced Design-Formula")
  
  # Verify covariates exist in metadata
  check_design_covariates(design_formula, metadata_df)
  
  # Make sure there are sufficient replicates for all factors
  check_replicates(metadata_df, all.vars(design_formula))
  
  # Add formulas to metadata for convenience
  attr(metadata_df, "design_formula") <- design_formula
  attr(metadata_df, "reduced_formula") <- reduced_formula
  
  return(metadata_df)
}

# Load and validate expression data
load_and_validate_expression_data <- function(args, metadata_df) {
  message("Loading expression data...")
  
  # Clean sample names for consistency
  clean_names <- clean_sample_names(args$name)
  
  # Load expression data
  expression_data_df <- load_expression_data(args$input, clean_names, READ_COL, RPKM_COL, INTERSECT_BY)
  message(glue::glue("Loaded expression data for {nrow(expression_data_df)} genes from {length(args$input)} files"))
  
  # Apply RPKM filtering if specified
  rpkm_filtered_count <- NULL
  if (!is.null(args$rpkm_cutoff)) {
    message(glue::glue("Applying RPKM cutoff of {args$rpkm_cutoff}..."))
    initial_gene_count <- nrow(expression_data_df)
    
    expression_data_df <- filter_rpkm(expression_data_df, args$rpkm_cutoff)
    
    # Ensure at least some genes remain
    if (nrow(expression_data_df) == 0) {
      report_error(
        "No genes remaining after RPKM filtering!",
        details = "All genes were removed due to the RPKM cutoff being too high.",
        recommendations = c("Reduce the RPKM threshold to retain more genes.")
      )
    }
    
    # Calculate how many genes were removed
    rpkm_filtered_count <- initial_gene_count - nrow(expression_data_df)
    
    message(glue::glue("{rpkm_filtered_count} genes removed, {nrow(expression_data_df)} genes retained"))
  }
  
  # Process count data
  read_counts_columns <- grep(
    paste(READ_COL, sep = ""),
    colnames(expression_data_df),
    value = TRUE,
    ignore.case = TRUE
  )
  
  read_counts_data_df <- expression_data_df %>%
    dplyr::mutate_at("GeneId", toupper) %>%
    dplyr::distinct(GeneId, .keep_all = TRUE) %>%
    remove_rownames() %>%
    column_to_rownames("GeneId")
  
  read_counts_data_df <- read_counts_data_df[read_counts_columns]
  
  # Clean count data column names
  colnames(read_counts_data_df) <- lapply(colnames(read_counts_data_df), function(s) {
    paste(head(unlist(strsplit(s, " ", fixed = TRUE)), -1), collapse = " ")
  })
  colnames(read_counts_data_df) <- clean_sample_names(colnames(read_counts_data_df))
  
  # Verify sample name consistency between metadata and counts
  validate_sample_consistency(metadata_df, read_counts_data_df)
  
  # Reorder count data columns to match metadata
  read_counts_data_df <- read_counts_data_df[, rownames(metadata_df)]
  
  # Return all results
  return(list(
    expression_df = expression_data_df,
    counts = read_counts_data_df,
    design_formula = attr(metadata_df, "design_formula"),
    reduced_formula = attr(metadata_df, "reduced_formula"),
    rpkm_filtered_count = rpkm_filtered_count
  ))
}

# Main execution function to organize workflow
run_deseq_analysis <- function(args) {
  # Start time tracking for the entire analysis
  total_start_time <- proc.time()
  
  # Setup and validation
  message("=== DESeq2 Analysis Pipeline ===")
  message(glue::glue("Analysis started at {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"))
  
  # 1. Load and validate metadata
  metadata_df <- load_and_validate_metadata(args)
  
  # 2. Load and validate expression data
  expression_data <- load_and_validate_expression_data(args, metadata_df)
  
  # 3. Process count data and apply batch correction
  count_data_results <- process_count_data(args, expression_data$counts, metadata_df, 
                                         expression_data$design_formula)
  
  # 4. Run DESeq2 analysis
  deseq_results <- run_deseq2(count_data_results$countData, 
                            metadata_df, 
                            count_data_results$design_formula,
                            args$reduced,
                            args)
  
  # 5. Process and export results
  export_results(deseq_results, expression_data$expression_df, 
                metadata_df, args, count_data_results$batch_warning,
                expression_data$rpkm_filtered_count)
  
  # Report total elapsed time
  total_elapsed <- proc.time() - total_start_time
  message(glue::glue("=== Analysis completed in {round(total_elapsed['elapsed']/60, 1)} minutes ==="))
  message(glue::glue("Finished at {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}"))
}

# Main script execution with enhanced error handling
main <- function() {
  # Set up error handling
  tryCatch({
    # Parse arguments
    args <- get_args()
    
    # Setup debugging and logging
    if ("debug" %in% names(args) && args$debug) {
      enable_debug_mode(args)
    }
    
    # Configure BiocParallel
    register(MulticoreParam(args$threads))
    log_message(glue::glue("Using {args$threads} CPU threads for parallel processing"), "INFO")
    
    # Run the analysis
    run_deseq_analysis(args)
    
    # Report successful completion
    log_message("DESeq2 analysis completed successfully! ✓", "SUCCESS")
    
  }, error = function(e) {
    # Capture full error information
    error_msg <- paste("Analysis failed:", e$message)
    error_trace <- paste(capture.output(traceback()), collapse="\n")
    
    log_message(error_msg, "ERROR")
    log_message("Error traceback:", "ERROR")
    message(error_trace)
    
    # Check if error_report.txt exists
    if (file.exists("error_report.txt")) {
      log_message("See error_report.txt for detailed error information", "ERROR")
    } else {
      # Create basic error report if none exists
      writeLines(
        paste0(
          "# ❌ Error Report\n\n",
          "**Timestamp:** ", format(Sys.time(), "%A, %B %d, %Y %I:%M:%S %p"), "\n\n",
          "**Error Message:** ", e$message, "\n\n",
          "## Error Traceback\n\n```\n", error_trace, "\n```\n"
        ),
        con = "error_report.txt"
      )
      log_message("Basic error report written to error_report.txt", "ERROR")
    }
    
    # Exit with error code
    quit(status = 1)
  })
} 