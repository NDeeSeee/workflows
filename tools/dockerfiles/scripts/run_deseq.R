#!/usr/bin/env Rscript
#
# Main entry point for DESeq Analysis
#

# Display startup message
message("Starting DESeq/DESeq2 Analysis")
message("Working directory:", getwd())

# Docker paths
workflow_file <- "/usr/local/bin/functions/deseq/workflow.R"
message("Looking for workflow file at:", workflow_file)

# Source the workflow file
if (file.exists(workflow_file)) {
  message("Found workflow file. Sourcing:", workflow_file)
  source(workflow_file)
  
  # Initialize the environment
  initialize_environment()
  
  # Run the workflow with memory management
  main_with_memory_management()
} else {
  message("ERROR: Workflow file not found at", workflow_file)
  message("Checking alternative locations...")
  
  # Try relative path
  relative_path <- "functions/deseq/workflow.R"
  if (file.exists(relative_path)) {
    message("Found workflow file at relative path:", relative_path)
    source(relative_path)
    initialize_environment()
    main_with_memory_management()
  } else {
    # Last resort - try to find it
    message("Attempting to locate workflow.R file...")
    system("find /usr/local -name workflow.R | grep deseq/", intern = FALSE)
    stop("Could not find workflow.R file. Please verify your installation.")
  }
}
