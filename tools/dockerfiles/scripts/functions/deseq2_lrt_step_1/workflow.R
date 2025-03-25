#!/usr/bin/env Rscript

# --- Main workflow functions ---

# Load and validate metadata
load_and_validate_metadata <- function(args) {
  message("Loading metadata...")
  
  # Get the file delimiter
  delimiter <- check_file_delimiter(args$meta)
  
  # Load metadata
  metadata_df <- read.table(
    args$meta,
    sep = delimiter,
    header = TRUE,
    stringsAsFactors = FALSE,
    row.names = 1
  )
  
  # Clean metadata column and row names
  colnames(metadata_df) <- clean_sample_names(colnames(metadata_df))
  rownames(metadata_df) <- clean_sample_names(rownames(metadata_df))
  
  message(glue::glue("Loaded metadata for {nrow(metadata_df)} samples with {ncol(metadata_df)} covariates"))
  
  # Check design formulas
  design_formula <- as.formula(args$design)
  
  # Apply comprehensive metadata validation using the common utility function
  metadata_df <- validate_metadata(metadata_df, args$batchcorrection, design_formula)
  
  # Add formulas to metadata for convenience
  attr(metadata_df, "design_formula") <- design_formula
  attr(metadata_df, "reduced_formula") <- as.formula(args$reduced)
  
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

# Harmonize parameter names for compatibility
harmonize_parameters <- function(args) {
  # Map renamed parameters to their original names for backwards compatibility
  param_mapping <- list(
    cluster_method = "cluster",
    row_distance = "rowdist", 
    column_distance = "columndist",
    k_hopach = "k",
    kmax_hopach = "kmax",
    output_prefix = "output"
  )
  
  # For each parameter, ensure we have a consistent name
  for (new_param in names(param_mapping)) {
    old_param <- param_mapping[[new_param]]
    
    # If the new parameter exists, use it
    if (!is.null(args[[new_param]])) {
      # Copy to old parameter name for functions that might use it
      args[[old_param]] <- args[[new_param]]
    } 
    # If only the old parameter exists, copy it to the new name
    else if (!is.null(args[[old_param]])) {
      args[[new_param]] <- args[[old_param]]
    }
  }
  
  return(args)
}

# Check required parameters for analysis
validate_analysis_params <- function(args) {
  # Check critical parameters
  critical_params <- c("input", "name", "meta", "design", "reduced")
  missing_params <- critical_params[!critical_params %in% names(args) | sapply(args[critical_params], is.null)]
  
  if (length(missing_params) > 0) {
    stop(paste("Missing critical parameters:", paste(missing_params, collapse=", ")))
  }
  
  # Validate batch correction parameter
  if (args$batchcorrection != "none" && !"batch" %in% colnames(read.table(args$meta, sep=get_file_type(args$meta), header=TRUE))) {
    warning("Batch correction requested but 'batch' column not found in metadata. Batch correction will be disabled.")
    args$batchcorrection <- "none"
  }
  
  return(args)
}

# Main wrapper function with memory management
main_with_memory_management <- function() {
  # Start timing
  start_time <- Sys.time()
  log_message("DESeq2 LRT Step 1 started", "START")
  
  # Get command line arguments
  args <- get_args()
  
  # Configure parallel processing
  if (args$threads > 1) {
    log_message(paste("Setting up parallel execution with", args$threads, "threads"), "CONFIG")
    register(MulticoreParam(args$threads))
  } else {
    log_message("Running in single-threaded mode", "CONFIG")
  }
  
  # Run the main workflow with validated args
  main(args)
  
  # Report end time and duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "secs")
  log_message(glue::glue("Total execution time: {round(as.numeric(duration), 2)} seconds"), "DONE")
  
  # Final memory report
  report_memory_usage("Final")
}

# Main script execution with enhanced error handling
main <- function(args = NULL) {
  # Set up error handling
  tryCatch({
    # Parse arguments if not provided
    if (is.null(args)) {
      args <- get_args()
    }
    
    # Harmonize parameter names for compatibility
    args <- harmonize_parameters(args)
    
    # Validate analysis parameters
    args <- validate_analysis_params(args)
    
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