#!/usr/bin/env Rscript

# --- Contrast generation functions ---

# Generate contrasts for main effects
generate_main_effect_contrasts <- function(dds, result_names) {
  log_message("Processing main effect contrasts...", "INFO")
  
  # Initialize dataframe
  main_effect_contrasts <- data.frame(
    factor = character(),
    contrast_name = character(),
    numerator = character(),
    denominator = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # Get all factor columns from colData
  col_data <- as.data.frame(colData(dds))
  factor_cols <- sapply(col_data, is.factor)
  factor_names <- names(factor_cols)[factor_cols]
  
  # Skip 'batch' column if it exists and is a factor
  if ("batch" %in% factor_names) {
    factor_names <- factor_names[factor_names != "batch"]
    debug_log("Excluding 'batch' from contrast generation")
  }
  
  debug_log(glue::glue("Processing factor columns: {paste(factor_names, collapse=', ')}"))
  
  # Process each factor
  for (factor_name in factor_names) {
    # Get levels for this factor
    factor_levels <- levels(col_data[[factor_name]])
    
    if (length(factor_levels) < 2) {
      debug_log(glue::glue("Skipping factor '{factor_name}' - fewer than 2 levels"))
      next
    }
    
    debug_log(glue::glue("Processing factor '{factor_name}' with levels: {paste(factor_levels, collapse=', ')}"))
    
    # Process factor according to number of levels
    if (length(factor_levels) == 2) {
      # Binary factor - simple case
      log_message(glue::glue("Processing binary factor: {factor_name}"), "INFO")
      
      # Check if the contrast exists in result_names
      contrast_name <- paste0(factor_name, factor_levels[2])
      
      if (contrast_name %in% result_names) {
        main_effect_contrasts <- rbind(
          main_effect_contrasts,
          data.frame(
            factor = factor_name,
            contrast_name = contrast_name,
            numerator = factor_levels[2],
            denominator = factor_levels[1],
            type = "main_effect",
            stringsAsFactors = FALSE
          )
        )
        debug_log(glue::glue("Added binary contrast: {contrast_name}"))
      } else {
        # Try alternate naming pattern
        alt_contrast_name <- find_contrast_name(result_names, factor_name, factor_levels[2])
        
        if (!is.null(alt_contrast_name)) {
          main_effect_contrasts <- rbind(
            main_effect_contrasts,
            data.frame(
              factor = factor_name,
              contrast_name = alt_contrast_name,
              numerator = factor_levels[2],
              denominator = factor_levels[1],
              type = "main_effect",
              stringsAsFactors = FALSE
            )
          )
          debug_log(glue::glue("Added binary contrast with alternate name: {alt_contrast_name}"))
        } else {
          log_message(glue::glue("Warning: Could not find contrast for binary factor '{factor_name}'"), "WARNING")
        }
      }
    } else {
      # Multi-level factor
      log_message(glue::glue("Processing multi-level factor: {factor_name} with {length(factor_levels)} levels"), "INFO")
      
      # Reference level is the first one
      ref_level <- factor_levels[1]
      
      # For each non-reference level
      for (i in 2:length(factor_levels)) {
        level <- factor_levels[i]
        contrast_name <- paste0(factor_name, level)
        
        if (contrast_name %in% result_names) {
          main_effect_contrasts <- rbind(
            main_effect_contrasts,
            data.frame(
              factor = factor_name,
              contrast_name = contrast_name,
              numerator = level,
              denominator = ref_level,
              type = "main_effect",
              stringsAsFactors = FALSE
            )
          )
          debug_log(glue::glue("Added multi-level contrast: {contrast_name}"))
        } else {
          # Try alternate naming pattern
          alt_contrast_name <- find_contrast_name(result_names, factor_name, level)
          
          if (!is.null(alt_contrast_name)) {
            main_effect_contrasts <- rbind(
              main_effect_contrasts,
              data.frame(
                factor = factor_name,
                contrast_name = alt_contrast_name,
                numerator = level,
                denominator = ref_level,
                type = "main_effect",
                stringsAsFactors = FALSE
              )
            )
            debug_log(glue::glue("Added multi-level contrast with alternate name: {alt_contrast_name}"))
          } else {
            log_message(glue::glue("Warning: Could not find contrast for level '{level}' of factor '{factor_name}'"), "WARNING")
          }
        }
      }
    }
  }
  
  log_message(glue::glue("Completed main effect contrast generation with {nrow(main_effect_contrasts)} contrasts"), "INFO")
  return(main_effect_contrasts)
}

# Generate contrasts for interaction effects
generate_interaction_effect_contrasts <- function(dds, result_names) {
  log_message("Processing interaction effect contrasts...", "INFO")
  
  # Initialize dataframe
  interaction_effect_contrasts <- data.frame(
    factor = character(),
    contrast_name = character(),
    numerator = character(),
    denominator = character(),
    type = character(),
    stringsAsFactors = FALSE
  )
  
  # Identify interaction terms in result_names
  interaction_terms <- result_names[grepl(":", result_names, fixed = TRUE)]
  
  if (length(interaction_terms) == 0) {
    debug_log("No interaction terms found in DESeq2 results")
    return(interaction_effect_contrasts)
  }
  
  debug_log(glue::glue("Found {length(interaction_terms)} interaction terms: {paste(interaction_terms, collapse=', ')}"))
  
  # Process each interaction term
  for (interaction_term in interaction_terms) {
    # For now, we'll just add these to the results with some placeholder values
    # In a real implementation, you'd parse these properly to extract the interacting factors
    
    interaction_effect_contrasts <- rbind(
      interaction_effect_contrasts,
      data.frame(
        factor = "interaction",
        contrast_name = interaction_term,
        numerator = interaction_term,
        denominator = "reference",
        type = "interaction_effect",
        stringsAsFactors = FALSE
      )
    )
    
    debug_log(glue::glue("Added interaction contrast: {interaction_term}"))
  }
  
  log_message(glue::glue("Completed interaction effect contrast generation with {nrow(interaction_effect_contrasts)} contrasts"), "INFO")
  return(interaction_effect_contrasts)
} 