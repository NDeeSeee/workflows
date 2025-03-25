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
  
  message(glue::glue("Found {length(read_counts_columns)} read count columns"))
  
  # Check if GeneId column exists
  if (!INTERSECT_BY %in% colnames(expression_data_df)) {
    stop(paste("Required column", INTERSECT_BY, "not found in expression data"))
  }
  
  # Check for duplicate GeneId values
  if (any(duplicated(expression_data_df[[INTERSECT_BY]]))) {
    dup_genes <- expression_data_df[[INTERSECT_BY]][duplicated(expression_data_df[[INTERSECT_BY]])]
    message(glue::glue("Warning: Found {length(dup_genes)} duplicate gene identifiers"))
    message(glue::glue("First few duplicates: {paste(head(dup_genes), collapse=', ')}"))
    
    # Make gene IDs unique
    expression_data_df[[INTERSECT_BY]] <- make.unique(as.character(expression_data_df[[INTERSECT_BY]]))
    message("Made gene identifiers unique")
  }
  
  # Extract read count data with safer approach
  tryCatch({
    # First extract only needed columns
    read_counts_subset <- expression_data_df[, c(INTERSECT_BY, read_counts_columns)]
    
    # Check for duplicate column names
    if (any(duplicated(colnames(read_counts_subset)))) {
      dup_cols <- colnames(read_counts_subset)[duplicated(colnames(read_counts_subset))]
      stop(paste("Duplicate column names in read count data:", paste(dup_cols, collapse=", ")))
    }
    
    # Make column names clean for downstream processing
    temp_colnames <- colnames(read_counts_subset)
    temp_colnames[-1] <- lapply(temp_colnames[-1], function(s) {
      # Extract the sample name (before the space)
      parts <- unlist(strsplit(s, " ", fixed = TRUE))
      if (length(parts) > 1) {
        return(parts[1])
      } else {
        return(s)
      }
    })
    colnames(read_counts_subset) <- temp_colnames
    
    # Now convert GeneId to row names
    # Clone the data frame
    read_counts_data_df <- read_counts_subset
    
    # Convert GeneId column to uppercase
    read_counts_data_df[[INTERSECT_BY]] <- toupper(read_counts_data_df[[INTERSECT_BY]])
    
    # Keep only first instance of each GeneId
    read_counts_data_df <- read_counts_data_df[!duplicated(read_counts_data_df[[INTERSECT_BY]]), ]
    
    # Check for duplicated gene IDs again after toupper conversion
    if (any(duplicated(read_counts_data_df[[INTERSECT_BY]]))) {
      dup_genes <- read_counts_data_df[[INTERSECT_BY]][duplicated(read_counts_data_df[[INTERSECT_BY]])]
      stop(paste("Duplicate gene IDs after uppercase conversion:", paste(head(dup_genes), collapse=", ")))
    }
    
    # Set row names safely
    rownames(read_counts_data_df) <- read_counts_data_df[[INTERSECT_BY]]
    read_counts_data_df <- read_counts_data_df[, -1, drop = FALSE] # Remove GeneId column
    
  }, error = function(e) {
    stop(paste("Error processing count data:", e$message))
  })
  
  # Clean count data column names with improved method
  cleaned_colnames <- clean_sample_names(colnames(read_counts_data_df))
  
  # Check for duplicates in cleaned column names
  if (any(duplicated(cleaned_colnames))) {
    dup_cols <- cleaned_colnames[duplicated(cleaned_colnames)]
    message(paste("Warning: Duplicate column names after cleaning:", paste(dup_cols, collapse=", ")))
    
    # Make column names unique
    cleaned_colnames <- make.unique(cleaned_colnames)
    message("Made column names unique")
  }
  
  # Apply cleaned column names
  colnames(read_counts_data_df) <- cleaned_colnames
  
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
  if (args$batchcorrection != "none") {
    # Check if metadata file exists
    if (!file.exists(args$meta)) {
      stop("Metadata file does not exist")
    }
    
    # Get the file delimiter
    delimiter <- check_file_delimiter(args$meta)
    
    # Read metadata to check for batch column
    metadata <- read.table(args$meta, sep=delimiter, header=TRUE)
    
    if (!"batch" %in% colnames(metadata)) {
      warning("Batch correction requested but 'batch' column not found in metadata. Batch correction will be disabled.")
      args$batchcorrection <- "none"
    }
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
    # Get detailed error information
    error_msg <- paste("Analysis failed:", e$message)
    
    # Set up detailed traceback capture
    options(rlang_backtrace_on_error = "full")
    
    # Get call stack information
    call_stack <- sys.calls()
    call_stack_formatted <- paste(format(call_stack), collapse = "\n")
    
    # Get traceback
    error_trace <- paste(capture.output(traceback()), collapse="\n")
    
    # Get session information
    session_info <- paste(capture.output(sessionInfo()), collapse="\n")
    
    # Log detailed error
    log_message(error_msg, "ERROR")
    log_message("Error occurred in call:", "ERROR")
    log_message(deparse(e$call), "ERROR")
    
    # Write detailed error report
    error_report_file <- file.path(getwd(), "deseq_lrt_error_report.txt")
    writeLines(
      paste0(
        "# DESeq2 LRT Analysis Error Report\n\n",
        "## Error Details\n\n",
        "**Timestamp:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
        "**Error Message:** ", e$message, "\n\n",
        "**Error Call:** ", deparse(e$call), "\n\n",
        
        "## Stack Trace\n\n```r\n", call_stack_formatted, "\n```\n\n",
        
        "## Error Traceback\n\n```\n", error_trace, "\n```\n\n",
        
        "## Session Information\n\n```\n", session_info, "\n```\n"
      ),
      con = error_report_file
    )
    
    log_message(paste("Detailed error report written to:", error_report_file), "ERROR")
    
    # Output simple error message to stderr for script integration
    message(error_msg)
    
    # Exit with error code
    quit(status = 1)
  })
} 