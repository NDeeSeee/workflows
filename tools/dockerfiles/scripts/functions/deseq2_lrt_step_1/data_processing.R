#!/usr/bin/env Rscript

# --- Data processing functions ---

# Function to scale data based on scaling_type
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