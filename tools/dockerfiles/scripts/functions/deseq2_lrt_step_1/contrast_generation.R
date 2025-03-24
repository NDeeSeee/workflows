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
  
  # Get factor information from colData
  col_data <- as.data.frame(colData(dds))
  factor_cols <- sapply(col_data, is.factor)
  factor_names <- names(factor_cols)[factor_cols]
  
  # Process each interaction term
  for (interaction_term in interaction_terms) {
    # Split the interaction term by colon
    interacting_factors <- unlist(strsplit(interaction_term, ":", fixed = TRUE))
    
    if (length(interacting_factors) < 2) {
      log_message(glue::glue("Invalid interaction term: {interaction_term}"), "WARNING")
      next
    }
    
    # Parse each factor and its level
    factor_levels <- list()
    reference_levels <- list()
    
    for (i in 1:length(interacting_factors)) {
      term <- interacting_factors[i]
      
      # Try different parsing approaches for factor_level
      if (grepl("_", term)) {
        # Format like "factor_level"
        parts <- unlist(strsplit(term, "_"))
        potential_factor <- parts[1]
        
        if (potential_factor %in% factor_names) {
          factor_levels[[i]] <- list(
            factor = potential_factor,
            level = paste(parts[-1], collapse = "_")
          )
          
          # Get reference level (first level)
          ref_levels <- levels(col_data[[potential_factor]])
          if (length(ref_levels) > 0) {
            reference_levels[[i]] <- ref_levels[1]
          } else {
            reference_levels[[i]] <- NA
          }
        }
      } else {
        # Try to match against known factors
        for (factor_name in factor_names) {
          if (startsWith(term, factor_name)) {
            level <- substring(term, nchar(factor_name) + 1)
            if (level != "") {
              factor_levels[[i]] <- list(
                factor = factor_name,
                level = level
              )
              
              # Get reference level
              ref_levels <- levels(col_data[[factor_name]])
              if (length(ref_levels) > 0) {
                reference_levels[[i]] <- ref_levels[1]
              } else {
                reference_levels[[i]] <- NA
              }
              break
            }
          }
        }
      }
      
      # If we couldn't determine the factor/level, log a warning
      if (is.null(factor_levels[[i]])) {
        log_message(glue::glue("Could not parse interaction term component: {term}"), "WARNING")
        factor_levels[[i]] <- list(factor = "unknown", level = term)
        reference_levels[[i]] <- NA
      }
    }
    
    # Build numerator and denominator strings
    num_factors <- sapply(factor_levels, function(x) x$factor)
    num_levels <- sapply(factor_levels, function(x) x$level)
    
    # Format as factor1level1:factor2level2
    numerator <- paste(
      paste0(num_factors, num_levels),
      collapse = ":"
    )
    
    # Format reference as factor1ref1:factor2ref2
    denominator <- paste(
      paste0(num_factors, reference_levels),
      collapse = ":"
    )
    
    # Add to results dataframe
    interaction_effect_contrasts <- rbind(
      interaction_effect_contrasts,
      data.frame(
        factor = paste(num_factors, collapse = ":"),
        contrast_name = interaction_term,
        numerator = numerator,
        denominator = denominator,
        type = "interaction_effect",
        stringsAsFactors = FALSE
      )
    )
    
    debug_log(glue::glue("Added interaction contrast: {interaction_term}"))
  }
  
  log_message(glue::glue("Completed interaction effect contrast generation with {nrow(interaction_effect_contrasts)} contrasts"), "INFO")
  return(interaction_effect_contrasts)
} 