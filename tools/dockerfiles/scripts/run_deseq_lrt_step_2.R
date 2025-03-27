#!/usr/bin/env Rscript
#
# Main entry point for DESeq2 LRT Step 2 Analysis
#

# Display startup message
message("Starting DESeq2 LRT Step 2 Analysis")
message("Working directory:", getwd())

# Docker paths
workflow_file <- "/usr/local/bin/functions/deseq2_lrt_step_2/workflow.R"
message("Looking for workflow file at:", workflow_file)

# Source the workflow file
if (file.exists(workflow_file)) {
  message("Found workflow file. Sourcing:", workflow_file)
  source(workflow_file)
  
  # Initialize the environment
  initialize_environment()
  
  # Get command line arguments
  args <- get_args()
  
  # Run the workflow
  run_workflow(args)
} else {
  message("ERROR: Workflow file not found at", workflow_file)
  message("Checking alternative locations...")
  
  # Try relative path
  relative_path <- "functions/deseq2_lrt_step_2/workflow.R"
  if (file.exists(relative_path)) {
    message("Found workflow file at relative path:", relative_path)
    source(relative_path)
    initialize_environment()
    args <- get_args()
    run_workflow(args)
  } else {
    # Last resort - try to find it
    message("Attempting to locate workflow.R file...")
    system("find /usr/local -name workflow.R | grep deseq2_lrt_step_2", intern = FALSE)
    stop("Could not find workflow.R file. Please verify your installation.")
  }
}
