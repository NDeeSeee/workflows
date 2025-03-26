#!/usr/bin/env Rscript

# --- Main workflow functions ---

# Load all required libraries for DESeq2 LRT analysis
load_required_libraries <- function() {
  suppressMessages({
    # Core packages
    library(argparse)
    library(BiocParallel)
    library(DESeq2)
    
    # Data manipulation
    library(tidyverse)
    library(data.table)
    library(conflicted)  # For resolving namespace conflicts
    
    # Batch correction
    library(limma)
    
    # Visualization
    library(pheatmap)
    library(RColorBrewer)
    library(ggplot2)
    library(ggrepel)
    library(plotly)
    
    # GCT export
    library(cmapR)
    
    # Utilities
    library(pryr)        # For memory usage tracking
    library(rlang)
    library(stringr)
    library(glue)
    library(logger)
  })
  
  message("All required libraries loaded successfully")
}

# Configure R options for DESeq2 analysis
configure_r_options <- function() {
  # Set warning level
  options(warn = -1)
  options(rlang_backtrace_on_error = "full")
  options("width" = 400)
  options(error = function() {
    message("An unexpected error occurred. Aborting script.")
    quit(save = "no", status = 1, runLast = FALSE)
  })

  # Set memory management options for large datasets
  options(future.globals.maxSize = 4000 * 1024^2)  # 4GB max for global data
  options(expressions = 5000)  # Increase expression stack size

  # Configure garbage collection behavior
  gcinfo(FALSE)  # Disable GC messages by default
  options(gc.aggressiveness = 0)  # Default GC behavior
  
  message("R options configured for DESeq2 analysis")
}

# Load all required source files
load_source_files <- function() {
  # Source utility functions from common directory
  # First try Docker standard path, then fall back to relative path
  if (file.exists("/usr/local/bin/functions/common/utilities.R")) {
    source("/usr/local/bin/functions/common/utilities.R")
  } else if (file.exists("functions/common/utilities.R")) {
    source("functions/common/utilities.R")
  } else {
    stop("Could not find utilities.R file")
  }

  report_memory_usage("Before loading function files")

  # Source common functions
  source_with_fallback("functions/common/visualization.R", "/usr/local/bin/functions/common/visualization.R")
  source_with_fallback("functions/common/export_functions.R", "/usr/local/bin/functions/common/export_functions.R")

  # Source DESeq2 LRT Step 1 specific functions
  source_with_fallback("functions/deseq2_lrt_step_1/cli_args.R", "/usr/local/bin/functions/deseq2_lrt_step_1/cli_args.R")
  source_with_fallback("functions/deseq2_lrt_step_1/data_processing.R", "/usr/local/bin/functions/deseq2_lrt_step_1/data_processing.R")
  source_with_fallback("functions/deseq2_lrt_step_1/deseq2_analysis.R", "/usr/local/bin/functions/deseq2_lrt_step_1/deseq2_analysis.R")
  source_with_fallback("functions/deseq2_lrt_step_1/contrast_generation.R", "/usr/local/bin/functions/deseq2_lrt_step_1/contrast_generation.R")

  report_memory_usage("After loading function files")
  
  # Configure plot theme
  configure_plot_theme()
  
  message("All source files loaded successfully")
}

# Pre-process command line arguments to handle unexpected positional arguments
preprocess_args <- function() {
  # Get all arguments
  all_args <- commandArgs(trailingOnly = TRUE)
  
  # Print original arguments for debugging
  message("Original arguments:")
  for (i in seq_along(all_args)) {
    message(paste("  ", i, ":", all_args[i]))
  }
  
  # Check if we have any arguments at all
  if (length(all_args) == 0) {
    message("No command line arguments found!")
    return(all_args)
  }
  
  # For array-based arguments, we need to preserve them during reconstruction
  array_args <- c("--input", "--name")
  array_values <- list()
  
  # Collect all array arguments first
  for (arg in array_args) {
    indices <- which(all_args == arg)
    if (length(indices) > 0) {
      values <- c()
      for (idx in indices) {
        if (idx < length(all_args) && !grepl("^--", all_args[idx+1])) {
          values <- c(values, all_args[idx+1])
        }
      }
      if (length(values) > 0) {
        array_values[[arg]] <- values
      }
    }
  }
  
  # Extract key required arguments first (meta, design, reduced)
  required_args <- c("--meta", "--design", "--reduced")
  required_values <- list()
  for (req_arg in required_args) {
    arg_index <- which(all_args == req_arg)
    if (length(arg_index) > 0 && arg_index[1] < length(all_args)) {
      required_values[[req_arg]] <- all_args[arg_index[1] + 1]
    }
  }
  
  # Extract all other singleton args (non-array args)
  other_args <- list()
  i <- 1
  while (i <= length(all_args)) {
    if (grepl("^--", all_args[i]) && !(all_args[i] %in% array_args)) {
      # If next arg exists and is not a flag
      if (i < length(all_args) && !grepl("^--", all_args[i+1])) {
        other_args[[all_args[i]]] <- all_args[i+1]
        i <- i + 2
      } else {
        # Just a flag
        other_args[[all_args[i]]] <- TRUE
        i <- i + 1
      }
    } else {
      # Skip array args (already handled) and positional args (handled later)
      i <- i + 1
    }
  }
  
  # Check if we have any arguments without -- prefix (positional arguments)
  positional_args <- all_args[!grepl("^--", all_args) & !grepl("^-[a-zA-Z]", all_args)]
  
  # Remove positional args that are values of flags
  flag_indices <- which(grepl("^--", all_args) | grepl("^-[a-zA-Z]", all_args))
  value_indices <- flag_indices + 1
  value_indices <- value_indices[value_indices <= length(all_args)]
  positional_values <- all_args[value_indices]
  positional_args <- setdiff(positional_args, positional_values)
  
  # If we have positional arguments, try to classify them
  if (length(positional_args) > 0) {
    message("Found positional arguments:", paste(positional_args, collapse=", "))
    
    # Separate positional arguments into names and inputs based on file path pattern
    file_pattern <- "\\.(tsv|csv)$"
    input_files <- positional_args[grepl(file_pattern, positional_args)]
    sample_names <- positional_args[!grepl(file_pattern, positional_args)]
    
    # If we have files or names, add them to our array_values collection
    if (length(input_files) > 0) {
      if ("--input" %in% names(array_values)) {
        array_values[["--input"]] <- c(array_values[["--input"]], input_files)
      } else {
        array_values[["--input"]] <- input_files
      }
    }
    
    if (length(sample_names) > 0) {
      if ("--name" %in% names(array_values)) {
        array_values[["--name"]] <- c(array_values[["--name"]], sample_names)
      } else {
        array_values[["--name"]] <- sample_names
      }
    }
    
    message("Reclassifying positional arguments:")
    message("  Input files:", paste(input_files, collapse=", "))
    message("  Sample names:", paste(sample_names, collapse=", "))
  }
  
  # Build new args list
  new_args <- c()
  
  # Add regular args
  for (arg_name in names(other_args)) {
    new_args <- c(new_args, arg_name, other_args[[arg_name]])
  }
  
  # Add required args if found
  for (req_arg in required_args) {
    if (req_arg %in% names(required_values)) {
      new_args <- c(new_args, req_arg, required_values[[req_arg]])
    }
  }
  
  # Add array args
  for (arg_name in names(array_values)) {
    for (value in array_values[[arg_name]]) {
      new_args <- c(new_args, arg_name, value)
    }
  }
  
  # Check if we have all the key arguments
  if (length(new_args) > 0) {
    message("Reformatted argument list:")
    i <- 1
    while (i <= length(new_args)) {
      if (i+1 <= length(new_args)) {
        message(paste("  ", new_args[i], new_args[i+1]))
        i <- i + 2
      } else {
        message(paste("  ", new_args[i], "TRUE"))
        i <- i + 1
      }
    }
    return(new_args)
  }
  
  # Return original arguments if no preprocessing needed
  return(all_args)
}

# Resolve namespace conflicts explicitly
resolve_namespace_conflicts <- function() {
  conflicted::conflict_prefer("filter", "dplyr")
  conflicted::conflict_prefer("select", "dplyr")
  conflicted::conflict_prefer("rename", "dplyr")
  conflicted::conflict_prefer("slice", "dplyr")
  conflicted::conflict_prefer("mutate", "dplyr")
  conflicted::conflict_prefer("group_by", "dplyr")
  conflicted::conflict_prefer("summarize", "dplyr")
  conflicted::conflict_prefer("arrange", "dplyr")
  conflicted::conflict_prefer("%>%", "magrittr")
  conflicted::conflict_prefer("%in%", "base")
  
  # Additional conflict resolutions
  conflicted::conflict_prefer("intersect", "base")
  conflicted::conflict_prefer("setdiff", "base")
  conflicted::conflict_prefer("union", "base")
  conflicted::conflict_prefer("as.data.frame", "base")
  conflicted::conflict_prefer("lag", "stats")
  conflicted::conflict_prefer("filter", "dplyr")
  conflicted::conflict_prefer("first", "dplyr")
  conflicted::conflict_prefer("last", "dplyr")
  conflicted::conflict_prefer("layout", "graphics")
  conflicted::conflict_prefer("plot", "graphics")
  conflicted::conflict_prefer("desc", "dplyr")
}

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
  # This function is retained for backward compatibility with any code
  # that might still expect the old parameter names
  # No action needed since we've handled parameter mapping in ArgumentParser
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

# Run the main process - adapter between main_with_memory_management and run_deseq_analysis
run_main_process <- function(args) {
  # Start timing
  start_time <- Sys.time()
  message("Starting DESeq2 LRT analysis process...")

  # Configure parallel processing based on thread count
  if (args$threads > 1) {
    log_message(paste("Setting up parallel execution with", args$threads, "threads"), "CONFIG")
    register(MulticoreParam(args$threads))
  } else {
    log_message("Running in single-threaded mode", "CONFIG")
  }

  # Run the main analysis workflow
  run_deseq_analysis(args)

  # Report completion and elapsed time
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  message(sprintf("DESeq2 LRT analysis completed in %.2f minutes", as.numeric(elapsed)))
}

# Initialize the DESeq2 analysis environment
initialize_environment <- function() {
  # Configure R options
  configure_r_options()
  
  # Load required libraries
  load_required_libraries()
  
  # Load source files
  load_source_files()
  
  # Resolve namespace conflicts
  resolve_namespace_conflicts()
  
  message("DESeq2 analysis environment initialized successfully")
}

# Main wrapper function with memory management
main_with_memory_management <- function() {
  # Start timing
  start_time <- Sys.time()
  log_message("DESeq2 LRT Step 1 started", "START")
  
  tryCatch({
    # Pre-process command line arguments first - improves compatibility with CWL
    processed_args <- preprocess_args()

    # Override command line arguments with processed ones
    # This is a hack to make the argparse parser work with our processed arguments
    if (!identical(processed_args, commandArgs(trailingOnly = TRUE))) {
      message("Command line arguments have been preprocessed for compatibility")
      
      # Create a safer version of the command line arguments
      safe_args <- c()
      for (arg in processed_args) {
        # Escape any quotes and special characters
        escaped_arg <- gsub('"', '\\\\"', arg)
        # Add quotes around arguments with spaces
        if (grepl(" ", escaped_arg)) {
          escaped_arg <- paste0('"', escaped_arg, '"')
        }
        safe_args <- c(safe_args, escaped_arg)
      }
      
      # Temporarily save the processed arguments to a file and source it
      temp_args_file <- tempfile(pattern = "args_", fileext = ".R")
      on.exit(unlink(temp_args_file), add = TRUE)
      
      cat("commandArgs_original <- commandArgs\n", file = temp_args_file)
      cat("commandArgs <- function(trailingOnly = FALSE) {\n", file = temp_args_file, append = TRUE)
      cat("  if (trailingOnly) {\n", file = temp_args_file, append = TRUE)
      
      # Convert the arguments to a proper R character vector
      arg_vector <- paste("c(", paste(shQuote(processed_args), collapse = ", "), ")")
      cat(paste0("    return(", arg_vector, ")\n"), file = temp_args_file, append = TRUE)
      
      cat("  } else {\n", file = temp_args_file, append = TRUE)
      cat("    return(commandArgs_original(FALSE))\n", file = temp_args_file, append = TRUE)
      cat("  }\n", file = temp_args_file, append = TRUE)
      cat("}\n", file = temp_args_file, append = TRUE)
      
      source(temp_args_file)
    }

    ### Parse CLI arguments
    message("Parsing command line arguments...")
    args <- tryCatch({
      params$get_cli_args()
    }, error = function(e) {
      message("Error during argument parsing: ", e$message)
      
      # Emergency fallback for required arguments
      message("EMERGENCY: Attempting to extract required arguments")
      
      # Get original command line args
      orig_args <- commandArgs(trailingOnly = TRUE)
      
      # Initialize empty result with NULL values
      emergency_args <- list(
        meta = NULL,
        design = NULL,
        reduced = NULL,
        input = character(0),
        name = character(0)
      )
      
      # Extract required arguments directly if parser failed
      # This is a more aggressive approach to handle malformed arguments
      required_flags <- c("--meta", "--design", "--reduced", "--input", "--name")
      
      # Handle regular args that aren't array-based
      for (flag in c("--meta", "--design", "--reduced")) {
        flag_idx <- which(orig_args == flag)
        if (length(flag_idx) > 0 && flag_idx[1] < length(orig_args)) {
          value <- orig_args[flag_idx[1] + 1]
          if (!grepl("^--", value)) {
            param_name <- gsub("^--", "", flag)
            emergency_args[[param_name]] <- value
            message("Found ", flag, " = ", value)
          }
        }
      }
      
      # Handle array args
      for (flag in c("--input", "--name")) {
        flag_indices <- which(orig_args == flag)
        values <- character(0)
        
        # Get all values for the flag
        for (idx in flag_indices) {
          if (idx < length(orig_args) && !grepl("^--", orig_args[idx + 1])) {
            values <- c(values, orig_args[idx + 1])
          }
        }
        
        # Store values if found
        if (length(values) > 0) {
          param_name <- gsub("^--", "", flag)
          emergency_args[[param_name]] <- values
          message("Found ", flag, " = ", paste(values, collapse=", "))
        }
      }
      
      # Check for positional arguments that could be inputs or sample names
      positional_args <- orig_args[!grepl("^--", orig_args)]
      positional_args <- setdiff(positional_args, unlist(emergency_args))
      
      if (length(positional_args) > 0) {
        message("Found ", length(positional_args), " potential positional arguments")
        
        # Try to identify input files vs sample names
        file_pattern <- "\\.(tsv|csv)$"
        input_files <- positional_args[grepl(file_pattern, positional_args)]
        sample_names <- positional_args[!grepl(file_pattern, positional_args)]
        
        if (length(input_files) > 0) {
          emergency_args$input <- c(emergency_args$input, input_files)
          message("Adding ", length(input_files), " input files from positional args")
        }
        
        if (length(sample_names) > 0) {
          emergency_args$name <- c(emergency_args$name, sample_names)
          message("Adding ", length(sample_names), " sample names from positional args")
        }
      }
      
      # Check if we found all required arguments
      missing_args <- names(emergency_args)[sapply(emergency_args, function(x) is.null(x) || (is.character(x) && length(x) == 0))]
      if (length(missing_args) > 0) {
        message("ERROR: Still missing required arguments: ", paste(missing_args, collapse=", "))
        stop("Missing required arguments: ", paste(missing_args, collapse=", "))
      }
      
      # Set defaults for other parameters
      default_params <- list(
        output = "results",
        rpkm_cutoff = 1,
        k_hopach = 5,
        kmax_hopach = 6,
        cosine_normalized = TRUE,
        distance_method = "cosine",
        clustering_method = "complete",
        filter = TRUE,
        batch_correct = FALSE,
        center_rowmeans = TRUE,
        fdr = 0.1
      )
      
      # Merge with existing args
      for (param in names(default_params)) {
        if (is.null(emergency_args[[param]])) {
          emergency_args[[param]] <- default_params[[param]]
        }
      }
      
      # Print final parameters
      message("Emergency argument parsing completed!")
      message("Parameters:")
      for (param in names(emergency_args)) {
        param_value <- emergency_args[[param]]
        if (is.character(param_value) && length(param_value) > 1) {
          message(paste0("  ", param, ": [", paste(param_value, collapse=", "), "]"))
        } else {
          message(paste0("  ", param, ": ", param_value))
        }
      }
      
      return(emergency_args)
    })
    
    ### Create output dir if not exist
    if (!dir.exists(args$output)) {
      dir.create(args$output, recursive = TRUE)
    }
    
    ### Setup logging
    logfile_path <- file.path(args$output, "log.txt")
    logger::log_appender(logger::appender_tee(logfile_path))
    set_log_level_verbose()
    
    log_info("Starting DESeq2 analysis with LRT")
    log_debug("Arguments:", print_all_args(args))
    
    ### Validate CLI arguments
    params$assert_args(args)
    log_debug("CLI arguments validated")
    
    ### Run main process
    run_main_process(args)
    
  }, error = function(e) {
    # Also log to stderr - important for debugging CWL
    write(paste("ERROR:", e$message), stderr())
    
    # Try to log to file if possible
    tryCatch({
      log_error(e$message)
      log_error(e$call)
      log_error(traceback())
    }, error = function(e) {
      write(paste("Failed to log error to file:", e$message), stderr())
    })
    
    # Make sure there's a non-zero exit 
    quit(status = 1)
  })
  
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