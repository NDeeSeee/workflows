#!/usr/bin/env Rscript
#
# Common utility functions for DESeq2 analysis
#

#' Configure a standard plot theme for consistent visualizations
#'
#' @export
configure_plot_theme <- function() {
  ggplot2::theme_set(
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.text = ggplot2::element_text(color = "black", size = 10),
      axis.title = ggplot2::element_text(color = "black", size = 12),
      legend.key = ggplot2::element_rect(fill = "white"),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black"),
      strip.text = ggplot2::element_text(color = "black", size = 10)
    )
  )
}

#' Min-max scaling function for data normalization
#'
#' @param x Numeric vector to scale
#' @param min_range Minimum value in output range
#' @param max_range Maximum value in output range
#' @return Scaled numeric vector
#' @export
scale_min_max <- function(x, min_range = -2, max_range = 2) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  
  if (max_val == min_val) {
    return(rep(0, length(x)))  # Avoid division by zero
  }
  
  scaled_x <- (x - min_val) / (max_val - min_val) * (max_range - min_range) + min_range
  return(scaled_x)
}

#' Filter expression data based on RPKM cutoff
#'
#' @param expression_df Data frame containing expression data
#' @param n RPKM cutoff value
#' @return Filtered data frame
#' @export
filter_rpkm <- function(expression_df, n) {
  expression_df %>%
    dplyr::filter(if_any(dplyr::contains("Rpkm"), ~. > n))
}

#' Extract factor names and levels from contrast string
#'
#' @param contrast_string String representation of the contrast
#' @return List with factor and level components
#' @export
extract_factors_and_levels <- function(contrast_string) {
  # Assuming format like "factor1_level1_vs_factor2_level2"
  parts <- strsplit(contrast_string, "_vs_")[[1]]
  
  if (length(parts) != 2) {
    stop("Contrast string does not match expected format: factor1_level1_vs_factor2_level2")
  }
  
  factor1_parts <- strsplit(parts[1], "_")[[1]]
  factor2_parts <- strsplit(parts[2], "_")[[1]]
  
  factor1 <- factor1_parts[1]
  factor2 <- factor2_parts[1]
  
  level1 <- paste(factor1_parts[-1], collapse = "_")
  level2 <- paste(factor2_parts[-1], collapse = "_")
  
  return(list(
    factor1 = factor1,
    level1 = level1,
    factor2 = factor2,
    level2 = level2
  ))
}

#' Validate metadata for DESeq2 analysis
#'
#' @param metadata_df Metadata data frame
#' @param batch_correction Batch correction method (or "none")
#' @param design_formula Design formula (optional)
#' @return Validated and processed metadata data frame
#' @export
validate_metadata <- function(metadata_df, batch_correction = "none", design_formula = NULL) {
  log_message("Validating metadata...", "STEP")
  
  # Check dimensions
  if (nrow(metadata_df) < 2) {
    report_error(
      "Insufficient samples in metadata",
      details = paste("Found only", nrow(metadata_df), "samples, minimum required is 2"),
      recommendations = c(
        "Ensure metadata file contains at least 2 samples"
      )
    )
  }
  
  # Check for empty or NA values
  for (col in colnames(metadata_df)) {
    na_count <- sum(is.na(metadata_df[[col]]))
    empty_count <- sum(metadata_df[[col]] == "", na.rm = TRUE)
    
    if (na_count > 0 || empty_count > 0) {
      report_error(
        "Missing values in metadata",
        details = paste0("Column '", col, "' contains ", 
                        na_count, " NA values and ", 
                        empty_count, " empty strings"),
        recommendations = c(
          "Fill in all missing values in the metadata file",
          "If certain values are truly unknown, use a placeholder value and note this in your analysis"
        )
      )
    }
  }
  
  # Convert character columns that should be factors to factors
  if (!is.null(design_formula)) {
    # Get variables from design formula
    design_vars <- all.vars(design_formula)
    
    # For each variable in the design formula
    for (var in design_vars) {
      # Skip if the variable doesn't exist
      if (!var %in% colnames(metadata_df)) {
        next
      }
      
      # If it's a character column, convert to factor
      if (is.character(metadata_df[[var]])) {
        log_message(paste0("Converting character column '", var, "' to factor"), "INFO")
        metadata_df[[var]] <- factor(metadata_df[[var]])
      }
    }
  }
  
  # Validate batch column if batch correction is requested
  if (batch_correction != "none") {
    # Check if batch column exists (case insensitive)
    batch_col <- colnames(metadata_df)[tolower(colnames(metadata_df)) == "batch"]
    
    if (length(batch_col) == 0) {
      report_error(
        "Batch column missing",
        details = "A column named 'batch' is required when batch correction is requested",
        recommendations = c(
          "Add a 'batch' column to your metadata file",
          "Or disable batch correction if not needed"
        )
      )
    }
    
    # Check batch column values
    if (!is.numeric(metadata_df[[batch_col]])) {
      # Try to convert to numeric if possible
      if (all(grepl("^[0-9]+$", metadata_df[[batch_col]]))) {
        log_message("Converting batch column to numeric", "WARNING")
        metadata_df[[batch_col]] <- as.numeric(metadata_df[[batch_col]])
      } else {
        report_error(
          "Invalid batch column values",
          details = "Batch column must contain numeric values",
          recommendations = c(
            "Encode batch information as numeric values (e.g., 1, 2, 3)",
            "Update metadata file with proper batch encoding"
          )
        )
      }
    }
    
    # Check that batch has at least 2 levels
    if (length(unique(metadata_df[[batch_col]])) < 2) {
      log_message("Warning: Only one batch detected. Batch correction may not be meaningful.", "WARNING")
    }
  }
  
  # Check for sufficient replicates in each group if design formula is provided
  if (!is.null(design_formula)) {
    validate_replicates(metadata_df, design_formula)
  }
  
  log_message("Metadata validation completed successfully", "SUCCESS")
  return(metadata_df)
}

#' Validate that there are sufficient replicates for DESeq2 analysis
#'
#' @param metadata_df Metadata data frame
#' @param design_formula Design formula
#' @return TRUE if validation passes, otherwise stops with error
#' @export
validate_replicates <- function(metadata_df, design_formula) {
  log_message("Checking for sufficient replicates...", "STEP")
  
  # Get all variables from the design formula
  formula_vars <- all.vars(design_formula)
  
  # Extract only the factors/grouping variables (skip batch)
  factor_vars <- formula_vars[formula_vars != "batch"]
  
  # For each factor variable
  for (var in factor_vars) {
    # Skip if not in metadata
    if (!var %in% colnames(metadata_df)) {
      log_message(paste0("Warning: Variable '", var, "' from design formula not found in metadata"), "WARNING")
      next
    }
    
    # Get the levels of this factor
    if (is.factor(metadata_df[[var]])) {
      var_levels <- levels(metadata_df[[var]])
    } else {
      var_levels <- unique(metadata_df[[var]])
    }
    
    # Check each level has at least 2 samples
    for (level in var_levels) {
      sample_count <- sum(metadata_df[[var]] == level)
      
      if (sample_count < 2) {
        report_error(
          "Insufficient replicates",
          details = paste0("Level '", level, "' of variable '", var, "' has only ", sample_count, " replicate(s)"),
          recommendations = c(
            "Each level/group should have at least 2 replicates for reliable differential expression analysis",
            "Consider removing or merging levels with insufficient replicates"
          )
        )
      }
    }
    
    log_message(paste0("Variable '", var, "' has sufficient replicates for all levels"), "SUCCESS")
  }
  
  # Check for interaction terms
  if (any(grepl(":", deparse(design_formula), fixed = TRUE))) {
    log_message("Checking replicates for interaction terms...", "INFO")
    
    # Extract interaction terms
    formula_str <- deparse(design_formula)
    interaction_terms <- unlist(stringr::str_extract_all(formula_str, "[^~+[:space:]]+:[^+[:space:]]+"))
    
    for (term in interaction_terms) {
      # Split the interaction term
      interacting_vars <- unlist(strsplit(term, ":", fixed = TRUE))
      
      # Check if both variables exist
      if (all(interacting_vars %in% colnames(metadata_df))) {
        # For each combination of levels
        var1_levels <- if (is.factor(metadata_df[[interacting_vars[1]]])) {
          levels(metadata_df[[interacting_vars[1]]])
        } else {
          unique(metadata_df[[interacting_vars[1]]])
        }
        
        var2_levels <- if (is.factor(metadata_df[[interacting_vars[2]]])) {
          levels(metadata_df[[interacting_vars[2]]])
        } else {
          unique(metadata_df[[interacting_vars[2]]])
        }
        
        # Check each interaction level has at least 2 samples
        for (level1 in var1_levels) {
          for (level2 in var2_levels) {
            sample_count <- sum(metadata_df[[interacting_vars[1]]] == level1 & 
                              metadata_df[[interacting_vars[2]]] == level2)
            
            if (sample_count < 2) {
              log_message(paste0("Warning: Interaction level '", level1, ":", level2, 
                                "' has only ", sample_count, " replicate(s)"), "WARNING")
            }
          }
        }
      }
    }
  }
  
  log_message("Replicate validation completed", "SUCCESS")
  return(TRUE)
}

#' Verify consistency between Step 1 metadata and Step 2 input
#'
#' @param step1_metadata Metadata from Step 1
#' @param step1_design_formula Design formula from Step 1
#' @param args Step 2 command-line arguments 
#' @return TRUE if consistent, otherwise stops with error
#' @export
verify_step_consistency <- function(step1_metadata, step1_design_formula, args) {
  log_message("Verifying consistency between Step 1 and Step 2...", "STEP")
  
  # Check if batch correction settings are compatible
  step1_batch_col <- colnames(step1_metadata)[tolower(colnames(step1_metadata)) == "batch"]
  has_batch_col <- length(step1_batch_col) > 0
  
  if (args$batchcorrection != "none" && !has_batch_col) {
    report_error(
      "Batch correction requested but no batch column found in Step 1 metadata",
      details = "Step 2 requested batch correction but Step 1 did not have a batch column",
      recommendations = c(
        "Set batchcorrection to 'none' in Step 2",
        "Or re-run Step 1 with batch information included in the metadata"
      )
    )
  }
  
  # Set batch warning if needed
  if (has_batch_col && args$batchcorrection == "none") {
    log_message("Warning: Batch column exists in metadata but batch correction is disabled", "WARNING")
  }
  
  # Verify all covariate levels from design formula exist in metadata
  if (!is.null(step1_design_formula)) {
    design_vars <- all.vars(step1_design_formula)
    
    for (var in design_vars) {
      if (var != "batch" && var %in% colnames(step1_metadata)) {
        if (!is.factor(step1_metadata[[var]])) {
          log_message(paste0("Warning: Design variable '", var, "' is not a factor in Step 1 metadata"), "WARNING")
        }
      }
    }
  }
  
  log_message("Step consistency verification completed", "SUCCESS")
  return(TRUE)
}

#' Apply ComBat-Seq batch correction to count data
#'
#' @param count_data_df Data frame of count data
#' @param metadata_df Data frame of sample metadata
#' @param batch_col Name of the batch column in metadata
#' @param design_formula Design formula for the model
#' @return Data frame of batch-corrected counts
#' @export
apply_combatseq_correction <- function(count_data_df, metadata_df, batch_col, design_formula) {
  log_message("Applying ComBat-Seq batch correction...", "STEP")
  
  # Extract batch information
  batch <- metadata_df[[batch_col]]
  
  # Get covariates from design formula (excluding batch)
  design_vars <- all.vars(design_formula)
  covar_cols <- design_vars[design_vars != "batch"] 
  
  # Only include covariates that are in metadata
  covar_cols <- intersect(covar_cols, colnames(metadata_df))
  
  # Create model matrix if covariates exist
  if (length(covar_cols) > 0) {
    # Create formula for covariates
    covar_formula <- as.formula(paste("~", paste(covar_cols, collapse = "+")))
    
    # Create model matrix
    covar_matrix <- model.matrix(covar_formula, data = metadata_df)[, -1, drop = FALSE]
    
    # Apply ComBat-Seq with covariates
    corrected_counts <- sva::ComBat_seq(
      counts = as.matrix(count_data_df),
      batch = batch,
      group = NULL,
      covar_mod = covar_matrix
    )
  } else {
    # Apply ComBat-Seq without covariates
    corrected_counts <- sva::ComBat_seq(
      counts = as.matrix(count_data_df),
      batch = batch
    )
  }
  
  # Convert back to data frame
  corrected_counts_df <- as.data.frame(corrected_counts)
  
  return(corrected_counts_df)
}

#' Apply limma-based batch correction to normalized expression data
#'
#' @param expression_data Normalized expression data matrix
#' @param metadata_df Data frame of sample metadata
#' @param batch_col Name of the batch column in metadata
#' @param design_matrix Design matrix for the model (optional)
#' @return Batch-corrected expression data matrix
#' @export
apply_limma_batch_correction <- function(expression_data, metadata_df, batch_col, design_matrix = NULL) {
  log_message("Applying limma-based batch correction to normalized data...", "STEP")
  
  # Extract batch information
  batch <- metadata_df[[batch_col]]
  
  # Create design matrix if not provided
  if (is.null(design_matrix)) {
    # Find factor columns to include in design
    factor_cols <- sapply(metadata_df, is.factor)
    factor_names <- names(factor_cols)[factor_cols & names(factor_cols) != batch_col]
    
    if (length(factor_names) > 0) {
      # Create formula for design
      design_formula <- as.formula(paste("~", paste(factor_names, collapse = "+")))
      design_matrix <- model.matrix(design_formula, data = metadata_df)
    } else {
      # If no factors available, use intercept-only model
      design_matrix <- model.matrix(~1, data = metadata_df)
    }
  }
  
  # Apply batch correction
  corrected_data <- limma::removeBatchEffect(
    expression_data,
    batch = batch,
    design = design_matrix
  )
  
  return(corrected_data)
}

#' Apply batch correction based on specified method
#'
#' @param count_data Count data matrix or data frame
#' @param metadata_df Data frame of sample metadata
#' @param batch_method Batch correction method ("none", "combatseq", or "limma")
#' @param design_formula Design formula (optional)
#' @param normalized Whether the data is already normalized (TRUE) or raw counts (FALSE)
#' @return Batch-corrected data
#' @export
apply_batch_correction <- function(count_data, metadata_df, batch_method = "none", 
                                  design_formula = NULL, normalized = FALSE) {
  if (batch_method == "none") {
    log_message("No batch correction applied", "INFO")
    return(count_data)
  }
  
  # Find batch column (accounting for case)
  batch_col <- colnames(metadata_df)[tolower(colnames(metadata_df)) == "batch"]
  
  if (length(batch_col) == 0) {
    log_message("No batch column found in metadata. Cannot apply batch correction.", "WARNING")
    return(count_data)
  }
  
  # Apply appropriate correction method
  if (batch_method == "combatseq" && !normalized) {
    # ComBat-Seq for raw count data
    corrected_data <- apply_combatseq_correction(count_data, metadata_df, batch_col, design_formula)
    log_message("ComBat-Seq batch correction applied to raw counts", "SUCCESS")
    
    # Apply rounding to counts (ComBat-Seq should produce integers, but round to be safe)
    corrected_data <- round(corrected_data)
    
  } else if (batch_method == "limma" || (batch_method == "combatseq" && normalized)) {
    # For normalized data or if limma specified, use limma's removeBatchEffect
    corrected_data <- apply_limma_batch_correction(count_data, metadata_df, batch_col)
    log_message("Limma batch correction applied to normalized data", "SUCCESS")
    
  } else {
    log_message(paste("Invalid batch correction configuration:", 
                     "method =", batch_method, 
                     "normalized =", normalized), "WARNING")
    return(count_data)
  }
  
  return(corrected_data)
}

# --- Utility functions ---

# Function to determine file type based on extension
get_file_type <- function(file_path) {
  # Get file extension
  ext <- tolower(tools::file_ext(file_path))
  
  # Return appropriate delimiter based on extension
  if (ext == "csv") {
    return(",")
  } else if (ext %in% c("tsv", "txt")) {
    return("\t")
  } else {
    # Default to tab if extension is unknown
    log_message(glue::glue("Warning: Unknown file extension '{ext}', defaulting to tab delimiter"), "WARNING")
    return("\t")
  }
}

# Function to check file delimiter
check_file_delimiter <- function(file_path) {
  # Get the first line of the file
  first_line <- readLines(file_path, n = 1)
  
  # Count occurrences of common delimiters
  comma_count <- stringr::str_count(first_line, ",")
  tab_count <- stringr::str_count(first_line, "\t")
  
  # Determine expected delimiter based on file extension
  expected_delimiter <- get_file_type(file_path)
  expected_name <- ifelse(expected_delimiter == ",", "comma", "tab")
  
  # Check if the delimiter usage matches the file extension
  if (expected_delimiter == "," && tab_count > comma_count) {
    report_error(
      "File format mismatch",
      details = c(
        "The file appears to be tab-delimited, but has a .csv extension.",
        paste0("First line contains ", tab_count, " tabs and ", comma_count, " commas.")
      ),
      recommendations = c(
        "Rename the file with a .tsv extension",
        "OR convert the file to use comma delimiters to match the .csv extension"
      )
    )
  } else if (expected_delimiter == "\t" && comma_count > tab_count) {
    report_error(
      "File format mismatch",
      details = c(
        "The file appears to be comma-delimited, but has a .tsv extension.",
        paste0("First line contains ", comma_count, " commas and ", tab_count, " tabs.")
      ),
      recommendations = c(
        "Rename the file with a .csv extension",
        "OR convert the file to use tab delimiters to match the .tsv extension"
      )
    )
  }
  
  # Warn if first row has no delimiters
  if (comma_count == 0 && tab_count == 0) {
    log_message(
      "Warning: First line of file contains no delimiters. This may indicate a formatting problem.",
      "WARNING"
    )
  }
}

# Function to clean sample names
clean_sample_names <- function(sample_names) {
  # Replace spaces and special characters with underscores
  cleaned_names <- gsub("[^a-zA-Z0-9]", "_", sample_names)
  
  # Remove duplicate underscores
  cleaned_names <- gsub("_+", "_", cleaned_names)
  
  # Remove leading/trailing underscores
  cleaned_names <- gsub("^_|_$", "", cleaned_names)
  
  return(cleaned_names)
}

# Utility function to remove row names from a data frame
remove_rownames <- function(df) {
  rownames(df) <- NULL
  return(df)
} 