#!/usr/bin/env Rscript

# --- Data processing functions ---

# Constants
READ_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")

# Function to load expression data
load_expression_data <- function(filenames,
                                 prefixes,
                                 read_colname,
                                 rpkm_colname,
                                 intersect_by) {
  collected_isoforms <- NULL
  for (i in 1:length(filenames)) {
    isoforms <- read.table(
      filenames[i],
      sep = get_file_type(filenames[i]),
      header = TRUE,
      stringsAsFactors = FALSE
    )
    # FOR TEST ONLY TO REDUCE NUMBER OF ROWS AND SPEED UP TESTING
    if (!is.null(args$test_mode) && args$test_mode) {
      print("Test mode is ON, each sample will be limited to 500 rows")
      isoforms <- isoforms %>%
        dplyr::slice_head(n = 500)
    } else {
      print("Test mode is OFF, processing all rows")
    }

    print(paste0("Load ", nrow(isoforms), " rows from ", filenames[i]))
    colnames(isoforms)[colnames(isoforms) == read_colname] <- paste(prefixes[i], read_colname, sep = " ")
    colnames(isoforms)[colnames(isoforms) == rpkm_colname] <- paste(prefixes[i], rpkm_colname, sep = " ")
    if (is.null(collected_isoforms)) {
      collected_isoforms <- isoforms
    } else {
      collected_isoforms <- dplyr::full_join(collected_isoforms, isoforms, by = intersect_by)
    }
  }
  print(paste0("Number of rows common for all loaded files ", nrow(collected_isoforms)))

  return(collected_isoforms)
}

# Function to filter by RPKM
filter_rpkm <- function(expression_data_df, rpkm_cutoff) {
  # Get all RPKM columns
  rpkm_columns <- grep(
    RPKM_COL,
    colnames(expression_data_df),
    value = TRUE,
    ignore.case = TRUE
  )
  
  # Ensure we have RPKM columns
  if (length(rpkm_columns) == 0) {
    report_error(
      "No RPKM columns found in expression data",
      details = c(
        "Expected column names containing '" + RPKM_COL + "' but none were found.",
        "Available columns:", 
        paste(colnames(expression_data_df), collapse = ", ")
      ),
      recommendations = c(
        "Ensure your expression data files include RPKM values.",
        "Check column naming in your input files."
      )
    )
  }
  
  # Calculate row means of RPKM
  rpkm_means <- rowMeans(expression_data_df[, rpkm_columns], na.rm = TRUE)
  
  # Filter genes based on RPKM cutoff
  filtered_data <- expression_data_df[rpkm_means >= rpkm_cutoff, ]
  
  return(filtered_data)
}

# Function to validate sample consistency
validate_sample_consistency <- function(metadata_df, count_data_df) {
  log_message("Validating sample consistency between metadata and expression data...", "STEP")
  
  debug_log("Metadata sample names:", rownames(metadata_df))
  debug_log("Count data column names:", colnames(count_data_df))
  
  # Check for exact matches
  samples_in_metadata_not_in_data <- setdiff(rownames(metadata_df), colnames(count_data_df))
  samples_in_data_not_in_metadata <- setdiff(colnames(count_data_df), rownames(metadata_df))
  
  # More advanced debugging - check for case-insensitive matches or partial matches
  if (length(samples_in_metadata_not_in_data) > 0) {
    debug_log("Checking for case-insensitive matches for missing samples...")
    
    # Try case-insensitive matching
    lower_count_names <- tolower(colnames(count_data_df))
    for (sample in samples_in_metadata_not_in_data) {
      matches <- grep(tolower(sample), lower_count_names, value = TRUE)
      if (length(matches) > 0) {
        debug_log(glue::glue("Sample '{sample}' has potential case-insensitive matches: {paste(matches, collapse=', ')}"))
      }
    }
  }
  
  # Report results
  if (length(samples_in_metadata_not_in_data) > 0) {
    report_error(
      "Missing Samples in Expression Data",
      details = c(
        "### Issue:",
        "The following samples are present in the metadata but missing from the expression data:",
        "",
        paste(paste0("- ", samples_in_metadata_not_in_data), collapse = "\n")
      ),
      recommendations = c(
        "Ensure that sample names in the metadata file match exactly those in the expression data.",
        "Check for typos, extra spaces, or differences in letter case.",
        "Verify that the correct expression data file was loaded."
      )
    )
  }
  
  if (length(samples_in_data_not_in_metadata) > 0) {
    report_error(
      "Missing Samples in Metadata",
      details = c(
        "### Issue:",
        "The following samples are present in the expression data but missing from the metadata:",
        "",
        paste(paste0("- ", samples_in_data_not_in_metadata), collapse = "\n")
      ),
      recommendations = c(
        "Ensure that all samples in the expression data are listed in the metadata file.",
        "Check for typos, extra spaces, or differences in letter case.",
        "Confirm that the metadata file corresponds to the current expression dataset."
      )
    )
  }
  
  log_message("Sample consistency validation passed", "SUCCESS")
} 