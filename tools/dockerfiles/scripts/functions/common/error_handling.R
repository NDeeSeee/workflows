#!/usr/bin/env Rscript

# --- Error handling functions ---

# Detailed error reporting function
report_error <- function(msg, details = NULL, recommendations = NULL) {
  # Format the timestamp in a friendly, readable way
  timestamp <- format(Sys.time(), "%A, %B %d, %Y %I:%M:%S %p")
  
  # Construct the Markdown report
  md_msg <- paste0(
    "# ❌ Error Report\n\n",
    "**Timestamp:** ", timestamp, "\n\n",
    "**Error Message:** ", msg, "\n\n"
  )
  
  # Add details section if provided
  if (!is.null(details)) {
    md_msg <- paste0(
      md_msg,
      "## Details\n\n",
      "```\n",
      paste(details, collapse = "\n"),
      "\n```\n\n"
    )
  }
  
  # Add recommendations section if provided
  if (!is.null(recommendations)) {
    recs <- paste0("- ", paste(recommendations, collapse = "\n- "))
    md_msg <- paste0(
      md_msg,
      "## Recommendations\n\n",
      recs, "\n"
    )
  }
  
  # Write the Markdown message to a file
  writeLines(md_msg, con = "error_report.txt")
  
  # Stop execution by raising an error with the given message
  stop(msg)
} 