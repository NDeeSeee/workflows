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
  
  # Source common utility functions
  source_with_fallback("functions/common/error_handling.R", "/usr/local/bin/functions/common/error_handling.R")
  source_with_fallback("functions/common/logging.R", "/usr/local/bin/functions/common/logging.R")
  source_with_fallback("functions/common/utilities.R", "/usr/local/bin/functions/common/utilities.R")
  
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

#' Helper function to source files with fallback paths
#'
#' @param filepath Relative path to file
#' @param absolute_path Absolute path to file (fallback)
#' @return Result of source() call
source_with_fallback <- function(filepath, absolute_path = NULL) {
  # Try absolute path first if provided
  if (!is.null(absolute_path) && file.exists(absolute_path)) {
    log_message(paste("Sourcing from absolute path:", absolute_path))
    return(source(absolute_path))
  }
  
  # Try relative path
  if (file.exists(filepath)) {
    log_message(paste("Sourcing from relative path:", filepath))
    return(source(filepath))
  }
  
  # Try a standard Docker path
  docker_path <- file.path("/usr/local/bin", filepath)
  if (file.exists(docker_path)) {
    log_message(paste("Sourcing from Docker path:", docker_path))
    return(source(docker_path))
  }
  
  # If all fails, error
  log_error(paste("Could not find file to source:", filepath))
  stop(paste("Could not find file to source:", filepath))
}

#' Report current memory usage
#'
#' @param label Label for the memory usage report
report_memory_usage <- function(label = "") {
  gc(verbose = FALSE)
  
  if (requireNamespace("pryr", quietly = TRUE)) {
    mem_used <- pryr::mem_used()
    log_message(paste0("[Memory] ", label, ": ", round(mem_used / 1024^2, 1), " MB"))
  } else {
    log_message(paste0("[Memory] ", label, ": Not available (pryr package not installed)"))
  }
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
    
    # Export MA plot
    export_ma_plot(
      deseq_results$res, 
      paste(args$output, "_ma_plot", sep = "")
    )
    
    # Export PCA plot
    export_pca_plot(
      deseq_results$vst, 
      paste(args$output, "_pca_plot", sep = ""),
      deseq_results$pca_intgroup
    )
    
    # Export MDS plot
    export_mds_html_plot(
      norm_counts_data = deseq_results$vst,
      location = paste(args$output, "mds_plot.html", sep = "_")
    )
    
    # Prepare matrix for heatmap
    vsd <- DESeq2::assay(deseq_results$vst)
    rownames(vsd) <- collected_isoforms[, c("GeneId")]
    mat <- get_top_expressed_genes(deseq_results$vst, deseq_results$norm_counts, column_data)
    
    # Process DESeq2 results
    DESeqRes <- process_deseq_results(
      deseq_results, 
      collected_isoforms, 
      read_count_cols, 
      args$digits
    )
    
    # Set normCounts for GCT export
    normCounts <- deseq_results$norm_counts
    rownames(normCounts) <- toupper(collected_isoforms[, c("GeneId")])
    
  } else {
    log_message("ERROR: Not enough replicates to run DESeq2 analysis")
    log_message("DESeq2 requires at least two replicates per condition")
    return(NULL)
  }
  
  report_memory_usage("After DESeq2 analysis")
  
  # Export expression heatmap
  export_heatmap(
    mat,
    column_data,
    paste(args$output, "_expression_heatmap", sep = "")
  )
  
  # Export DESeq results to TSV file
  results_filename <- paste(args$output, "_report.tsv", sep = "")
  write.table(
    DESeqRes,
    file = results_filename,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  log_message(paste("Exported DESeq report to", results_filename))
  
  # Process data for GCT export
  log_message("Preparing data for GCT export")
  
  # Prepare row metadata for GCT
  row_metadata <- DESeqRes %>%
    dplyr::mutate_at("GeneId", toupper) %>%
    dplyr::distinct(GeneId, .keep_all = TRUE) %>%
    dplyr::remove_rownames() %>%
    dplyr::column_to_rownames("GeneId") %>%
    dplyr::select(log2FoldChange, pvalue, padj) %>%
    dplyr::arrange(desc(log2FoldChange))
  
  # Prepare column metadata for GCT
  col_metadata <- column_data %>%
    dplyr::mutate_at(colnames(.), as.vector)
  
  # Export normalized counts to GCT format
  export_gct(
    counts_mat = normCounts,
    row_metadata = row_metadata,
    col_metadata = col_metadata,
    location = paste(args$output, "_counts_all.gct", sep = "")
  )
  
  # Get size of matrix before filtering
  read_count_matrix_all_size <- dim(normCounts)
  log_message(paste("Size of normalized counts matrix before filtering:", 
                    read_count_matrix_all_size[1], "rows by", 
                    read_count_matrix_all_size[2], "columns"))
  
  # Filter by significance
  log_message(paste("Filtering normalized counts to include only features with padj <=", args$fdr))
  row_metadata <- row_metadata %>%
    dplyr::filter(.$padj <= args$fdr)
  
  # Apply log2FoldChange threshold if needed
  if (args$lfcthreshold > 0) {
    log_message(paste("Applying log2FoldChange threshold of", args$lfcthreshold))
    
    if (args$regulation == "up") {
      row_metadata <- row_metadata %>%
        dplyr::filter(.$log2FoldChange >= args$lfcthreshold)
    } else if (args$regulation == "down") {
      row_metadata <- row_metadata %>%
        dplyr::filter(.$log2FoldChange <= -args$lfcthreshold)
    } else {
      row_metadata <- row_metadata %>%
        dplyr::filter(abs(.$log2FoldChange) >= args$lfcthreshold)
    }
  }
  
  # Subset normCounts based on filtering
  normCounts <- normCounts[as.vector(rownames(row_metadata)), , drop = FALSE]
  log_message(paste("Size of normalized counts matrix after filtering:", 
                    nrow(normCounts), "rows by", 
                    ncol(normCounts), "columns"))
  
  # Perform clustering if requested
  if (!is.null(args$cluster)) {
    log_message(paste("Performing", args$cluster, "clustering"))
    
    k <- args$k
    kmax <- args$kmax
    
    # Column clustering if requested
    if (args$cluster == "column" || args$cluster == "both") {
      log_message("Clustering columns")
      clustered_data_cols <- get_clustered_data(
        normCounts, 
        by = "col", 
        k = k, 
        kmax = kmax, 
        scaling_type = args$scaling_type, 
        dist = args$columndist
      )
      normCounts <- normCounts[, clustered_data_cols$order, drop = FALSE]
      col_metadata <- col_metadata[clustered_data_cols$order, , drop = FALSE]
      col_metadata <- cbind(col_metadata, clustered_data_cols$clusters)
    }
    
    # Row clustering if requested
    if (args$cluster == "row" || args$cluster == "both") {
      log_message("Clustering rows")
      clustered_data_rows <- get_clustered_data(
        normCounts, 
        by = "row", 
        k = k, 
        kmax = kmax, 
        scaling_type = args$scaling_type, 
        dist = args$rowdist
      )
      normCounts <- clustered_data_rows$expression[clustered_data_rows$order, , drop = FALSE]
      row_metadata <- row_metadata[clustered_data_rows$order, , drop = FALSE]
      row_metadata <- cbind(row_metadata, clustered_data_rows$clusters)
    }
  }
  
  # Export filtered normalized counts to GCT format
  log_message("Exporting filtered normalized counts to GCT format")
  export_gct(
    counts_mat = normCounts,
    row_metadata = row_metadata,
    col_metadata = col_metadata,
    location = paste(args$output, "_counts_filtered.gct", sep = "")
  )
  
  # Export CLS phenotype data
  log_message("Exporting CLS phenotype data")
  export_cls(
    categories = col_metadata[, "conditions"], 
    paste(args$output, "_phenotypes.cls", sep = "")
  )
  
  report_memory_usage("End of workflow")
  log_message("DESeq workflow completed successfully")
  
  # Return results
  return(list(
    DESeqRes = DESeqRes,
    normCounts = normCounts,
    row_metadata = row_metadata,
    col_metadata = col_metadata
  ))
}

#' Main function with memory management
#'
#' @return Result of workflow execution
main_with_memory_management <- function() {
  # Start timing
  start_time <- Sys.time()
  log_message(sprintf("DESeq analysis started at %s", format(start_time, "%Y-%m-%d %H:%M:%S")))
  
  # Get command line arguments
  args <- with_error_handling({
    get_args()
  })
  
  if (is.null(args)) {
    log_error("Failed to parse command line arguments")
    stop("Failed to parse command line arguments")
  }
  
  # Configure parallel processing
  if (args$threads > 1) {
    log_message(paste("Setting up parallel execution with", args$threads, "threads"))
    BiocParallel::register(BiocParallel::MulticoreParam(args$threads))
  } else {
    log_message("Running in single-threaded mode")
  }
  
  # Run the main workflow
  results <- with_error_handling({
    run_workflow(args)
  })
  
  if (is.null(results)) {
    log_error("Workflow execution failed")
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  # Report end time and duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  log_message(sprintf("DESeq analysis completed at %s", format(end_time, "%Y-%m-%d %H:%M:%S")))
  log_message(sprintf("Total execution time: %.2f minutes", round(as.numeric(duration), 2)))
  
  # Clean up large objects to free memory
  rm(results)
  invisible(gc())
  
  # Final memory report
  report_memory_usage("Final")
  
  return(TRUE)
} 