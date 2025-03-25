#!/usr/bin/env Rscript

# --- Data processing functions ---

# Constants for column name patterns
READ_COL <- "Read"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- "GeneId"

# Function to load expression data from multiple files
load_expression_data <- function(input_files, sample_names, read_col="Read", rpkm_col="Rpkm", intersect_by="GeneId") {
  # Process each input file
  expression_data_list <- list()
  
  for (i in 1:length(input_files)) {
    input_file <- input_files[i]
    sample_name <- sample_names[i]
    
    # Determine file type
    delimiter <- check_file_delimiter(input_file)
    
    # Load data
    message(paste("Loading file:", input_file))
    
    data <- read.table(
      input_file, 
      sep = delimiter,
      header = TRUE, 
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      comment.char = "#"
    )
    
    # Rename columns to include sample name
    data_cols <- colnames(data)
    for (col in data_cols) {
      if (grepl(read_col, col, ignore.case = TRUE) || grepl(rpkm_col, col, ignore.case = TRUE)) {
        # Rename to include sample name
        new_col <- paste(sample_name, col, sep = " ")
        colnames(data)[colnames(data) == col] <- new_col
      }
    }
    
    expression_data_list[[i]] <- data
  }
  
  # Merge data by gene ID
  merged_data <- Reduce(function(x, y) {
    merge(x, y, by = intersect_by, all = TRUE, sort = FALSE)
  }, expression_data_list)
  
  # Replace NA values with 0
  merged_data[is.na(merged_data)] <- 0
  
  return(merged_data)
}

# Function to filter expression data based on RPKM values
filter_rpkm <- function(expression_data, rpkm_cutoff) {
  if (is.null(rpkm_cutoff) || rpkm_cutoff <= 0) {
    message("No RPKM filtering applied")
    return(expression_data)
  }
  
  message(paste("Filtering data with RPKM cutoff:", rpkm_cutoff))
  
  # Get RPKM columns
  rpkm_cols <- grep(RPKM_COL, colnames(expression_data), ignore.case = TRUE, value = TRUE)
  
  if (length(rpkm_cols) == 0) {
    warning("No RPKM columns found, skipping RPKM filtering")
    return(expression_data)
  }
  
  # Keep rows where any RPKM value exceeds cutoff
  keep_rows <- apply(expression_data[, rpkm_cols, drop = FALSE], 1, function(x) {
    any(x > rpkm_cutoff, na.rm = TRUE)
  })
  
  filtered_data <- expression_data[keep_rows, , drop = FALSE]
  
  message(paste("RPKM filtering: removed", sum(!keep_rows), "rows, kept", sum(keep_rows), "rows"))
  
  return(filtered_data)
}

# Function to process count data and apply batch correction
process_count_data <- function(args, counts, metadata, design_formula) {
  message("Processing count data...")
  
  # Set up variables for batch correction result
  batch_warning <- NULL
  countData <- counts
  
  # Apply batch correction if requested
  if (args$batchcorrection != "none") {
    message(paste("Applying batch correction method:", args$batchcorrection))
    
    if (!"batch" %in% colnames(metadata)) {
      batch_warning <- "Batch correction requested but no 'batch' column found in metadata. Batch correction skipped."
      warning(batch_warning)
    } else {
      # Check batch correction method
      if (args$batchcorrection == "combatseq") {
        # Apply ComBat-Seq batch correction
        countData <- apply_combatseq_correction(counts, metadata, design_formula)
      } else if (args$batchcorrection == "model") {
        # Include batch in model (no data modification needed)
        design_formula <- update(design_formula, ~ . + batch)
        message("Added batch to design formula: ", deparse(design_formula))
      } else {
        batch_warning <- paste("Unknown batch correction method:", args$batchcorrection, ". Batch correction skipped.")
        warning(batch_warning)
      }
    }
  }
  
  # Return processed data
  return(list(
    countData = countData,
    design_formula = design_formula,
    batch_warning = batch_warning
  ))
}

# Helper function for ComBat-Seq batch correction
apply_combatseq_correction <- function(counts, metadata, design_formula) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    warning("Package 'sva' is required for ComBat-Seq batch correction but not available. Skipping batch correction.")
    return(counts)
  }
  
  message("Applying ComBat-Seq batch correction...")
  
  # Extract batch information
  batch <- metadata$batch
  
  # Get model matrix from design formula
  mod <- model.matrix(design_formula, data = metadata)
  
  # Apply ComBat-Seq
  corrected_counts <- sva::ComBat_seq(
    counts = as.matrix(counts),
    batch = batch,
    group = mod
  )
  
  # Convert back to data frame
  corrected_counts <- as.data.frame(corrected_counts)
  
  # Ensure row and column names are preserved
  rownames(corrected_counts) <- rownames(counts)
  colnames(corrected_counts) <- colnames(counts)
  
  message("Batch correction applied successfully")
  
  return(corrected_counts)
}

# Function to scale expression data based on scaling_type
scale_expression_data <- function(expression_data, scaling_type) {
  if (scaling_type == "none" || is.null(scaling_type)) {
    message("No scaling applied to expression data")
    return(expression_data)
  }
  
  message(paste("Applying", scaling_type, "scaling to expression data"))
  
  # Apply scaling
  if (scaling_type == "minmax") {
    # Min-max scaling to range [-2, 2]
    scaled_data <- apply(expression_data, 1, function(x) {
      min_val <- min(x, na.rm = TRUE)
      max_val <- max(x, na.rm = TRUE)
      if (max_val == min_val) {
        return(rep(0, length(x)))
      }
      return(4 * (x - min_val) / (max_val - min_val) - 2)
    })
    # Transpose back to original orientation
    scaled_data <- t(scaled_data)
    
  } else if (scaling_type == "zscore") {
    # Z-score scaling (mean=0, sd=1)
    scaled_data <- t(apply(expression_data, 1, function(x) {
      mean_val <- mean(x, na.rm = TRUE)
      sd_val <- sd(x, na.rm = TRUE)
      if (sd_val == 0) {
        return(rep(0, length(x)))
      }
      return((x - mean_val) / sd_val)
    }))
  } else {
    stop(paste("Unsupported scaling type:", scaling_type))
  }
  
  return(scaled_data)
}

# Function to perform clustering based on cluster method
perform_clustering <- function(expression_data, cluster_method, row_distance, column_distance, k, kmax) {
  if (cluster_method == "none" || is.null(cluster_method)) {
    message("No clustering performed")
    return(list(expression_data = expression_data, row_order = 1:nrow(expression_data), col_order = 1:ncol(expression_data)))
  }
  
  message(paste("Performing", cluster_method, "clustering"))
  
  # Load required library
  require(hopach)
  
  # Initialize clustering results
  row_clustering <- NULL
  col_clustering <- NULL
  
  # Perform row clustering if requested
  if (cluster_method %in% c("row", "both")) {
    message(paste("Row clustering with distance:", row_distance))
    row_clustering <- hopach::hopach(
      data = expression_data,
      dmat = NULL,
      dist.method = row_distance,
      k = k,
      kmax = kmax,
      khigh = NULL
    )
  }
  
  # Perform column clustering if requested
  if (cluster_method %in% c("column", "both")) {
    message(paste("Column clustering with distance:", column_distance))
    col_clustering <- hopach::hopach(
      data = t(expression_data),  # Transpose for column clustering
      dmat = NULL,
      dist.method = column_distance,
      k = k,
      kmax = kmax,
      khigh = NULL
    )
  }
  
  # Reorder data based on clustering
  row_order <- if (!is.null(row_clustering)) row_clustering$final$order else 1:nrow(expression_data)
  col_order <- if (!is.null(col_clustering)) col_clustering$final$order else 1:ncol(expression_data)
  
  ordered_data <- expression_data[row_order, col_order]
  
  return(list(
    expression_data = ordered_data,
    row_order = row_order,
    col_order = col_order,
    row_clustering = row_clustering,
    col_clustering = col_clustering
  ))
}

# Function to filter expression data based on RPKM cutoff
filter_by_rpkm <- function(expression_data, sample_metadata, rpkm_cutoff) {
  if (is.null(rpkm_cutoff)) {
    message("No RPKM-based filtering applied")
    return(expression_data)
  }
  
  message(paste("Filtering data with RPKM cutoff:", rpkm_cutoff))
  
  # Find columns containing "Rpkm" in their names
  rpkm_cols <- grep("Rpkm", colnames(expression_data), ignore.case = TRUE)
  
  if (length(rpkm_cols) == 0) {
    warning("No columns with 'Rpkm' in their names found, no filtering applied")
    return(expression_data)
  }
  
  # Keep rows where any RPKM column exceeds the cutoff
  keep_rows <- apply(expression_data[, rpkm_cols, drop = FALSE], 1, function(x) {
    any(x > rpkm_cutoff, na.rm = TRUE)
  })
  
  message(paste("Filtered out", sum(!keep_rows), "of", nrow(expression_data), "rows"))
  
  return(expression_data[keep_rows, ])
} 