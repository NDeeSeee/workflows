#!/usr/bin/env Rscript
#
# Main workflow for DESeq/DESeq2 differential expression analysis
#
# This file orchestrates the entire workflow for DESeq analysis,
# from environment setup to results processing and export.
#
# Version: 0.1.0

#' Initialize the environment for DESeq analysis
#'
#' This function loads required libraries, sources dependency files,
#' and sets up error handling and logging.
initialize_environment <- function() {
  # Set options
  options(warn = -1)
  options(rlang_backtrace_on_error = "full")
  options("width" = 400)
  options(error = function() {
    message("An unexpected error occurred. Aborting script.")
    quit(save = "no", status = 1, runLast = FALSE)
  })
  
  # Set memory management options for large datasets
  options(future.globals.maxSize = 4000 * 1024^2)  # 4GB max for global data
  options(expressions = 5000)  # Increase expression stack size
  
  # Configure garbage collection behavior
  gcinfo(FALSE)  # Disable GC messages by default
  options(gc.aggressiveness = 0)  # Default GC behavior
  
  # First load utilities which has source_with_fallback
  if (file.exists("functions/common/utilities.R")) {
    source("functions/common/utilities.R")
  } else if (file.exists("/usr/local/bin/functions/common/utilities.R")) {
    source("/usr/local/bin/functions/common/utilities.R")
  } else {
    stop("Could not find common utilities.R file, which is required")
  }
  
  # Source common utility functions
  source_with_fallback("functions/common/constants.R", "/usr/local/bin/functions/common/constants.R")
  source_with_fallback("functions/common/error_handling.R", "/usr/local/bin/functions/common/error_handling.R")
  source_with_fallback("functions/common/logging.R", "/usr/local/bin/functions/common/logging.R")
  
  # Source common visualization and export functions
  source_with_fallback("functions/common/visualization.R", "/usr/local/bin/functions/common/visualization.R")
  source_with_fallback("functions/common/clustering.R", "/usr/local/bin/functions/common/clustering.R")
  source_with_fallback("functions/common/export_functions.R", "/usr/local/bin/functions/common/export_functions.R")
  
  # Source DESeq-specific functions
  source_with_fallback("functions/deseq/cli_args.R", "/usr/local/bin/functions/deseq/cli_args.R")
  source_with_fallback("functions/deseq/data_processing.R", "/usr/local/bin/functions/deseq/data_processing.R")
  source_with_fallback("functions/deseq/deseq_analysis.R", "/usr/local/bin/functions/deseq/deseq_analysis.R")
  
  # Load required libraries
  load_required_libraries()
  
  # Configure plot theme
  configure_plot_theme()
  
  # Log initialization
  log_message("Environment initialized for DESeq analysis")
}

#' Load all required libraries for DESeq analysis
load_required_libraries <- function() {
  log_message("Loading required libraries")
  
  suppressMessages({
    # For argument parsing
    library(argparse)
    
    # For parallel processing
    library(BiocParallel)
    
    # For DESeq2 analysis
    library(DESeq2)
    
    # For data manipulation
    library(tidyverse)
    library(data.table)
    
    # For batch correction
    library(limma)
    
    # For clustering
    library(hopach)
    
    # For visualization
    library(pheatmap)
    library(RColorBrewer)
    library(gridExtra)
    library(ggplot2)
    library(ggrepel)
    library(Glimma)
    
    # For GCT export
    library(cmapR)
    
    # For memory profiling
    if (requireNamespace("pryr", quietly = TRUE)) {
      library(pryr)
    }
  })
  
  # Define dplyr functions with proper namespace to avoid conflicts
  `%>%` <- magrittr::`%>%`
  `%in%` <- base::`%in%`
  `%/%` <- base::`%/%`
  `%%` <- base::`%%`
  
  log_message("Libraries loaded successfully")
}

#' Configure plot theme for consistent visualization
configure_plot_theme <- function() {
  log_message("Configuring plot theme")
  
  # Set default theme for ggplot2
  theme_set(theme_bw() + 
            theme(text = element_text(size = 12),
                  axis.text = element_text(size = 10),
                  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  legend.title = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  strip.background = element_rect(fill = "lightgray"),
                  strip.text = element_text(face = "bold")))
}

#' Main workflow function
#'
#' @param args Command line arguments
#' @return Results of the analysis
run_workflow <- function(args) {
  log_message("Starting DESeq workflow")
  
  # Load isoforms/genes/tss files
  report_memory_usage("Before loading data")
  
  raw_data <- load_isoform_set(
    args$treated,
    args$talias,
    READ_COL,
    RPKM_COL,
    RPKM_TREATED_ALIAS,
    args$tname,
    INTERSECT_BY,
    args$digits,
    args$batchfile,
    load_isoform_set(
      args$untreated,
      args$ualias,
      READ_COL,
      RPKM_COL,
      RPKM_UNTREATED_ALIAS,
      args$uname,
      INTERSECT_BY,
      args$digits,
      args$batchfile
    )
  )
  
  # Extract data components
  collected_isoforms <- raw_data$collected_isoforms
  read_count_cols <- raw_data$read_colnames
  column_data <- raw_data$column_data
  
  log_message(paste("Number of rows common for all input files:", nrow(collected_isoforms)))
  
  # Apply RPKM filtering if requested
  if (!is.null(args$rpkm_cutoff)) {
    collected_isoforms <- filter_rpkm(collected_isoforms, args$rpkm_cutoff)
    log_message(paste("Expression data after RPKM filtering:", nrow(collected_isoforms), "rows"))
  }
  
  # Prepare count data for DESeq2
  count_data <- collected_isoforms[read_count_cols]
  log_message("Count data prepared for DESeq2 analysis")
  
  report_memory_usage("After loading data")
  
  # Run DESeq or DESeq2 based on sample count
  if (length(args$treated) > 1 && length(args$untreated) > 1) {
    log_message("Running DESeq2 analysis (multiple replicates available)")
    
    # Define design formula
    if (!is.null(args$batchfile) && args$batchcorrection == "model") {
      design <- ~conditions + batch
      log_message("Using design formula with batch effect: ~conditions + batch")
      batch_data <- args$batchfile$batch
    } else {
      design <- ~conditions
      log_message("Using standard design formula: ~conditions")
      batch_data <- NULL
    }
    
    # Run DESeq2 analysis
    deseq_results <- run_deseq2_analysis(
      count_data = count_data,
      col_data = column_data,
      design = design,
      batch_correction = args$batchcorrection,
      batch_data = batch_data,
      condition_names = c(condition1 = args$uname, condition2 = args$tname),
      args = args
    )
    
    # Generate summary markdown
    generate_md(
      args$batchcorrection, 
      args$batchfile, 
      deseq_results$res, 
      paste0(args$output, "_summary.md")
    )
    
    # Create summary plots
    create_summary_plots(
      deseq_results$dds,
      deseq_results$res,
      args$output,
      args$vst,
      args$pval,
      args$lfc
    )
    
    # Export data in various formats
    export_data(
      deseq_results$res,
      deseq_results$norm_counts,
      collected_isoforms,
      args
    )
    
    # Generate additional visualizations
    generate_visualizations(deseq_results, args)
    
    report_memory_usage("After DESeq2 analysis")
    
    # Return results
    return(deseq_results)
    
  } else {
    log_message("Running DESeq analysis with EdgeR (single replicate mode)")
    
    # Run EdgeR-based analysis for single replicate mode
    deseq_results <- run_edger_analysis(
      count_data = count_data,
      col_data = column_data,
      condition_names = c(condition1 = args$uname, condition2 = args$tname),
      args = args
    )
    
    # Export data in various formats
    export_data(
      deseq_results$res,
      deseq_results$norm_counts,
      collected_isoforms,
      args
    )
    
    # Return results
    return(deseq_results)
  }
}

#' Generate and save visualizations from DESeq2 results
#'
#' @param deseq_results Results from DESeq2 analysis
#' @param args Command line arguments
#' @return None
generate_visualizations <- function(deseq_results, args) {
  log_message("Generating visualizations")
  
  # Extract components from results
  dds <- deseq_results$dds
  res <- deseq_results$res
  
  # Create MA plot
  pdf(paste0(args$output, "_ma_plot.pdf"), width = 8, height = 6)
  DESeq2::plotMA(res, main = "MA Plot", ylim = c(-5, 5))
  dev.off()
  
  # Create dispersion plot
  pdf(paste0(args$output, "_dispersion_plot.pdf"), width = 8, height = 6)
  DESeq2::plotDispEsts(dds, main = "Dispersion Estimates")
  dev.off()
  
  # Create PCA plot if transformed data is available
  if (!is.null(deseq_results$vst_data)) {
    vst_data <- deseq_results$vst_data
    
    pdf(paste0(args$output, "_pca_plot.pdf"), width = 10, height = 8)
    DESeq2::plotPCA(vst_data, intgroup = "conditions") +
      ggtitle("PCA Plot") +
      theme(plot.title = element_text(hjust = 0.5))
    dev.off()
  }
  
  log_message("Visualizations saved successfully")
}

#' Wrapper function with memory management
#'
#' @return Result of workflow execution
main_with_memory_management <- function() {
  # Start timing
  start_time <- Sys.time()
  log_message(paste("DESeq analysis started at", format(start_time, "%Y-%m-%d %H:%M:%S")))
  
  # Get command line arguments
  args <- get_args()
  
  # Configure parallel processing if requested
  if (!is.null(args$threads) && args$threads > 1) {
    log_message(paste("Setting up parallel execution with", args$threads, "threads"))
    register(MulticoreParam(args$threads))
  } else {
    log_message("Running in single-threaded mode")
  }
  
  # Run the main workflow
  results <- run_workflow(args)
  
  # Report end time and duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  log_message(paste("DESeq analysis completed at", format(end_time, "%Y-%m-%d %H:%M:%S")))
  log_message(paste("Total execution time:", round(as.numeric(duration), 2), "minutes"))
  
  # Clean up large objects to free memory
  rm(results)
  invisible(gc())
  
  # Final memory report
  report_memory_usage("Final")
  
  log_message("DESeq analysis completed successfully")
} 