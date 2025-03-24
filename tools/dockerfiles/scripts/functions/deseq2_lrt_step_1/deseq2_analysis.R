#!/usr/bin/env Rscript

# --- DESeq2 analysis functions ---

# Process count data and apply batch correction
process_count_data <- function(args, count_data_df, metadata_df, design_formula) {
  log_message("Processing count data...", "STEP")
  
  # Initialize batch warning flag
  batch_warning <- FALSE
  
  # Check if batch correction is requested and apply it if needed
  if (!is.null(args$batchcorrection) && args$batchcorrection != "none") {
    log_message(glue::glue("Applying {args$batchcorrection} batch correction..."), "STEP")
    
    # Apply batch correction using common utility
    corrected_count_data <- apply_batch_correction(
      count_data = count_data_df,
      metadata_df = metadata_df,
      batch_method = args$batchcorrection,
      design_formula = design_formula,
      normalized = FALSE  # Raw counts for DESeq2
    )
    
    # Check if batch correction was actually applied
    if (identical(corrected_count_data, count_data_df)) {
      batch_warning <- TRUE
      log_message("Batch correction could not be applied to raw counts. Will be applied after normalization.", "WARNING")
    } else {
      log_message(paste("Batch correction applied using", args$batchcorrection), "SUCCESS")
    }
  } else {
    # No batch correction
    corrected_count_data <- count_data_df
    log_message("No batch correction applied", "INFO")
  }
  
  # Return processed data
  return(list(
    countData = corrected_count_data,
    design_formula = design_formula,
    batch_warning = batch_warning
  ))
}

# Run DESeq2 LRT analysis
run_deseq2 <- function(count_data, metadata_df, design_formula, reduced_formula, args) {
  log_message("Running DESeq2 LRT analysis...", "STEP")
  
  # Create DESeq2 dataset object
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = metadata_df,
    design = design_formula
  )
  
  # Pre-filter low count genes to speed up analysis
  log_message(glue::glue("Pre-filtering genes with fewer than {args$mincounts} counts across all samples..."), "INFO")
  keep <- rowSums(counts(dds)) >= args$mincounts
  dds <- dds[keep, ]
  log_message(glue::glue("Retained {sum(keep)} out of {length(keep)} genes after pre-filtering"), "INFO")
  
  # Store original design for reference
  design_formula_str <- deparse(design_formula)
  reduced_formula_str <- deparse(as.formula(args$reduced))
  
  # Verify formulas are different
  if (design_formula_str == reduced_formula_str) {
    report_error(
      "Full and reduced design formulas are identical",
      details = c(
        "Full design formula: ", design_formula_str,
        "Reduced design formula: ", reduced_formula_str
      ),
      recommendations = c(
        "The full design formula should include terms that are not in the reduced formula",
        "Modify either the full or reduced formula to create a proper comparison"
      )
    )
  }
  
  # Run DESeq with LRT test
  log_message("Running DESeq2 with Likelihood Ratio Test...", "STEP")
  dds <- DESeq2::DESeq(
    dds,
    test = "LRT",
    reduced = as.formula(args$reduced),
    parallel = TRUE
  )
  
  # Get LRT results
  lrt_res <- DESeq2::results(dds)
  
  # Generate contrasts
  contrasts <- generate_contrasts(dds, args)
  
  # Get normalized counts
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  
  log_message("DESeq2 analysis completed successfully", "SUCCESS")
  
  # Return results
  return(list(
    dds = dds,
    lrt_res = lrt_res,
    contrasts = contrasts,
    normCounts = norm_counts
  ))
}

# Generate all possible contrasts from DESeq2 results object
generate_contrasts <- function(dds, args) {
  log_message("Generating contrasts for differential expression analysis...", "STEP")
  
  # Get results names from DESeq2 object
  result_names <- resultsNames(dds)
  debug_log(glue::glue("Available result names from DESeq2: {paste(result_names, collapse=', ')}"))
  
  # Skip intercept
  result_names <- result_names[result_names != "Intercept"]
  
  # Initialize contrast dataframe  
  contrast_df <- data.frame(
    factor = character(),
    contrast_name = character(),
    numerator = character(),
    denominator = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # If no valid results, return empty dataframe
  if (length(result_names) == 0) {
    log_message("No valid contrasts found in DESeq2 results.", "WARNING")
    return(contrast_df)
  }
  
  # Process main effects
  main_effect_contrasts <- generate_main_effect_contrasts(dds, result_names)
  if (nrow(main_effect_contrasts) > 0) {
    contrast_df <- rbind(contrast_df, main_effect_contrasts)
    log_message(glue::glue("Generated {nrow(main_effect_contrasts)} main effect contrasts"), "SUCCESS")
  } else {
    log_message("No main effect contrasts generated", "INFO")
  }
  
  # Process interaction effects
  interaction_effect_contrasts <- generate_interaction_effect_contrasts(dds, result_names)
  if (nrow(interaction_effect_contrasts) > 0) {
    contrast_df <- rbind(contrast_df, interaction_effect_contrasts)
    log_message(glue::glue("Generated {nrow(interaction_effect_contrasts)} interaction effect contrasts"), "SUCCESS")
  } else {
    log_message("No interaction effect contrasts generated", "INFO")
  }
  
  # Final validation
  if (nrow(contrast_df) == 0) {
    log_message("Warning: No contrasts could be generated. Check your design formula and metadata.", "WARNING")
  } else {
    log_message(glue::glue("Successfully generated {nrow(contrast_df)} total contrasts"), "SUCCESS")
  }
  
  return(contrast_df)
}

# Helper function to find contrast name using various patterns
find_contrast_name <- function(result_names, factor_name, level) {
  # Common patterns for contrast names in DESeq2
  patterns <- c(
    paste0(factor_name, "_", level),
    paste0(factor_name, level),
    paste0(factor_name, ".", level),
    paste0(tolower(factor_name), "_", tolower(level)),
    paste0(tolower(factor_name), tolower(level)),
    paste0(tolower(factor_name), ".", tolower(level))
  )
  
  debug_log(glue::glue("Searching for contrast patterns: {paste(patterns, collapse=', ')}"))
  
  # Try exact matches first
  for (pattern in patterns) {
    if (pattern %in% result_names) {
      return(pattern)
    }
  }
  
  # Try partial/fuzzy matches
  for (pattern in patterns) {
    matches <- grep(pattern, result_names, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      debug_log(glue::glue("Found fuzzy match: {matches[1]} for pattern {pattern}"))
      return(matches[1])
    }
  }
  
  # No match found
  return(NULL)
} 