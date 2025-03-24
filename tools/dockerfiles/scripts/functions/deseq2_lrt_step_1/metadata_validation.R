#!/usr/bin/env Rscript

# --- Metadata validation functions ---

# Check design formula
check_design_formula <- function(formula_str, formula_name) {
  log_message(glue::glue("Validating {formula_name}..."), "STEP")
  
  tryCatch(
    {
      formula_obj <- as.formula(tolower(formula_str))
      log_message(glue::glue("{formula_name} is valid: {deparse(formula_obj)}"), "SUCCESS")
      return(formula_obj)
    },
    error = function(e) {
      report_error(
        glue::glue("Invalid {formula_name}: {formula_str}"),
        details = c(
          "The formula string could not be parsed as a valid R formula.",
          paste("Error:", e$message)
        ),
        recommendations = c(
          "Ensure the formula follows R syntax: ~factor1 + factor2",
          "Check for typos and matching parentheses",
          "Variables in formula should match column names in metadata",
          "Examples of valid formulas: ~condition, ~condition + batch, ~condition*time"
        )
      )
    }
  )
}

# Check covariates in design formula exist in metadata
check_design_covariates <- function(formula_obj, metadata_df) {
  log_message("Validating covariates in design formula...", "STEP")
  
  # Get all variables from the formula
  covariates <- all.vars(formula_obj)
  
  # Check if all covariates exist in metadata columns
  missing_covariates <- setdiff(covariates, colnames(metadata_df))
  
  if (length(missing_covariates) > 0) {
    # Get close matches for each missing covariate to help with debugging
    close_matches <- list()
    for (cov in missing_covariates) {
      # Find column names that partially match
      potential_matches <- grep(cov, colnames(metadata_df), value = TRUE, ignore.case = TRUE)
      if (length(potential_matches) > 0) {
        close_matches[[cov]] <- potential_matches
      }
    }
    
    # Prepare error details with suggestions
    details <- c(
      "The following covariates in your design formula are missing from metadata:",
      paste("- ", missing_covariates, collapse = "\n")
    )
    
    if (length(close_matches) > 0) {
      details <- c(
        details,
        "",
        "Potential similar column names in your metadata:"
      )
      
      for (cov in names(close_matches)) {
        details <- c(
          details,
          paste0("For '", cov, "', did you mean: ", paste(close_matches[[cov]], collapse = ", "))
        )
      }
    }
    
    # Add available columns to help debugging
    details <- c(
      details,
      "",
      "Available columns in metadata:",
      paste(colnames(metadata_df), collapse = ", ")
    )
    
    report_error(
      "Missing covariates in metadata",
      details = details,
      recommendations = c(
        "Ensure all variables in your design formula match column names in your metadata",
        "Check for typos or case sensitivity issues",
        "Add the missing columns to your metadata if they are required"
      )
    )
  }
  
  log_message("All covariates in design formula exist in metadata", "SUCCESS")
}

# Check batch column if batch correction is requested
check_batch_column <- function(metadata_df, batch_correction) {
  if (!is.null(batch_correction) && batch_correction != "none") {
    log_message("Checking batch column for batch correction...", "STEP")
    
    if (!"batch" %in% tolower(colnames(metadata_df))) {
      report_error(
        "Missing 'batch' column in metadata",
        details = c(
          "Batch correction was requested, but no 'batch' column was found in metadata.",
          "",
          "Available columns in metadata:",
          paste(colnames(metadata_df), collapse = ", ")
        ),
        recommendations = c(
          "Add a 'batch' column to your metadata file",
          "OR use a different column for batch information by renaming it to 'batch'",
          "OR disable batch correction by setting --batchcorrection to 'none'"
        )
      )
    }
    
    # Find the actual column name (accounting for case differences)
    batch_col <- colnames(metadata_df)[tolower(colnames(metadata_df)) == "batch"]
    
    # Check if batch column has at least two levels
    batch_values <- unique(metadata_df[[batch_col]])
    
    if (length(batch_values) < 2) {
      report_error(
        "Insufficient batch diversity",
        details = c(
          paste("The batch column contains only", length(batch_values), "unique value(s):"),
          paste(batch_values, collapse = ", ")
        ),
        recommendations = c(
          "Ensure the batch column contains at least two different batch values",
          "OR disable batch correction by setting --batchcorrection to 'none'"
        )
      )
    }
    
    log_message(glue::glue("Batch column validated with {length(batch_values)} unique batches"), "SUCCESS")
  }
}

# Check if there are sufficient replicates for all factors
check_replicates <- function(metadata_df, factor_names) {
  log_message("Checking for sufficient biological replicates...", "STEP")
  
  insufficient_replicates <- FALSE
  error_details <- c("The following conditions have insufficient replicates (at least 2 required):")
  
  for (factor_name in factor_names) {
    # Skip if factor doesn't exist in metadata (this should be caught by other checks)
    if (!factor_name %in% colnames(metadata_df)) {
      next
    }
    
    # Convert to factor if it's not already
    if (!is.factor(metadata_df[[factor_name]])) {
      metadata_df[[factor_name]] <- as.factor(metadata_df[[factor_name]])
    }
    
    # Count replicates for each level
    factor_levels <- levels(metadata_df[[factor_name]])
    
    for (level in factor_levels) {
      # Count samples with this factor level
      replicate_count <- sum(metadata_df[[factor_name]] == level)
      
      if (replicate_count < 2) {
        insufficient_replicates <- TRUE
        error_details <- c(
          error_details,
          paste0("- Factor '", factor_name, "', level '", level, "': only ", replicate_count, " sample(s)")
        )
      }
    }
  }
  
  if (insufficient_replicates) {
    report_error(
      "Insufficient biological replicates",
      details = error_details,
      recommendations = c(
        "DESeq2 requires at least 2 biological replicates per condition for reliable results",
        "Either add more samples to your dataset",
        "OR remove the under-replicated factor levels from your analysis",
        "OR consider using a different analysis method for low-replicate data"
      )
    )
  }
  
  log_message("All factors have sufficient replicates", "SUCCESS")
} 