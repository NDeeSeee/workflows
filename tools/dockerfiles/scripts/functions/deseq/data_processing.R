#!/usr/bin/env Rscript
#
# Data processing functions for DESeq/DESeq2 differential expression analysis
#
# This file contains functions for loading, preprocessing, and filtering
# expression data for the DESeq analysis workflow.
#
# Version: 0.1.0

# Define constants for column names
READ_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")
RPKM_UNTREATED_ALIAS <- "RpkmCondition1"
RPKM_TREATED_ALIAS <- "RpkmCondition2"

#' Determine file type from filename extension
#'
#' @param filename The path to the file
#' @return The appropriate separator character for the file type
get_file_type <- function(filename) {
  ext <- tools::file_ext(filename)
  separator <- ","
  if (ext == "tsv") {
    separator <- "\t"
  }
  return(separator)
}

#' Filter expression data by RPKM threshold
#'
#' @param expression_df Data frame containing expression data
#' @param n RPKM threshold value
#' @return Filtered data frame
filter_rpkm <- function(expression_df, n) {
  log_message(paste("Filtering expression data with RPKM cutoff:", n))
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for RPKM filtering")
  }
  
  filtered_df <- expression_df %>%
    dplyr::filter(dplyr::if_any(dplyr::contains("Rpkm"), ~. > n))
  
  log_message(paste("Filtered", nrow(expression_df) - nrow(filtered_df), "rows by RPKM cutoff"))
  
  return(filtered_df)
}

#' Load expression data from isoform files and merge them
#'
#' @param filenames Vector of file paths
#' @param prefixes Vector of sample prefixes
#' @param read_colname Column name for read counts
#' @param rpkm_colname Column name for RPKM values
#' @param rpkm_colname_alias Alias for combined RPKM values
#' @param conditions Condition name (e.g., "treated" or "untreated")
#' @param intersect_by Vector of column names to use for merging
#' @param digits Number of digits for rounding
#' @param batch_metadata Batch metadata for multi-factor analysis
#' @param collected_data Previously collected data (optional)
#' @return List containing merged expression data and metadata
load_isoform_set <- function(filenames,
                             prefixes,
                             read_colname,
                             rpkm_colname,
                             rpkm_colname_alias,
                             conditions,
                             intersect_by,
                             digits,
                             batch_metadata,
                             collected_data = NULL) {
  
  # Initialize rpkm_colnames if needed
  if (is.null(collected_data)) {
    rpkm_colnames <- character(0)
  } else {
    rpkm_colnames <- collected_data$rpkm_colnames
  }
  
  # Process each input file
  for (i in 1:length(filenames)) {
    # Load data
    isoforms <- with_error_handling({
      read.table(
        filenames[i],
        sep = get_file_type(filenames[i]),
        header = TRUE,
        stringsAsFactors = FALSE
      )
    })
    
    if (is.null(isoforms)) {
      log_error(paste("Failed to load file:", filenames[i]))
      next
    }
    
    # Rename columns
    new_read_colname <- paste(prefixes[i], conditions, sep = "_")
    colnames(isoforms)[colnames(isoforms) == read_colname] <- new_read_colname
    colnames(isoforms)[colnames(isoforms) == rpkm_colname] <- paste(conditions, i, rpkm_colname, sep = " ")
    
    # Create column data
    if (!is.null(batch_metadata)) {
      batch <- batch_metadata[prefixes[i], "batch"]
      log_message(
        paste0(
          "Loaded ", nrow(isoforms),
          " rows from '", filenames[i],
          "' as '", new_read_colname,
          "', batch '", batch, "'"
        )
      )
      column_data_frame <- data.frame(conditions, batch, row.names = c(new_read_colname))
    } else {
      log_message(
        paste0(
          "Loaded ", nrow(isoforms),
          " rows from '", filenames[i],
          "' as '", new_read_colname, "'"
        )
      )
      column_data_frame <- data.frame(conditions, row.names = c(new_read_colname))
    }
    
    # Initialize or update collected data
    if (is.null(collected_data)) {
      collected_data <- list(
        collected_isoforms = isoforms,
        read_colnames = c(new_read_colname),
        column_data = column_data_frame,
        rpkm_colnames = character(0)
      )
    } else {
      collected_data$collected_isoforms <- merge(
        collected_data$collected_isoforms,
        isoforms,
        by = intersect_by,
        sort = FALSE
      )
      collected_data$read_colnames <- c(collected_data$read_colnames, new_read_colname)
      collected_data$column_data <- rbind(collected_data$column_data, column_data_frame)
    }
  }
  
  # Calculate average RPKM for the condition
  rpkm_columns <- grep(
    paste("^", conditions, " [0-9]+ ", rpkm_colname, sep = ""),
    colnames(collected_data$collected_isoforms),
    value = TRUE,
    ignore.case = TRUE
  )
  
  collected_data$collected_isoforms[rpkm_colname_alias] <- format(
    rowSums(collected_data$collected_isoforms[, rpkm_columns, drop = FALSE]) / length(filenames),
    digits = digits
  )
  
  # Update RPKM column names and remove individual RPKM columns
  collected_data$rpkm_colnames <- c(collected_data$rpkm_colnames, rpkm_colname_alias)
  collected_data$collected_isoforms <- collected_data$collected_isoforms[, !colnames(collected_data$collected_isoforms) %in% rpkm_columns]
  
  return(collected_data)
}

#' Generate an analysis summary markdown file
#'
#' @param batchcorrection The batch correction method used
#' @param batchfile The batch file provided (if any)
#' @param deseq_results Results from DESeq/DESeq2 analysis
#' @param output_file Path to the output markdown file
#' @return None
generate_md <- function(batchcorrection, batchfile, deseq_results, output_file) {
  # Initialize the markdown content
  md_content <- ""

  # Add warning message if applicable
  if (!is.null(batchcorrection) && (batchcorrection == "combatseq" || batchcorrection == "model") && is.null(batchfile)) {
    warning_message <- "# Warning!\n\n---\n\n**You provided a batch-correction method, but not a batch-file.**\n\nThe chosen parameter was ignored.\n\nPlease ensure that you provide a batch file when using the following batch correction methods:\n\n- **combatseq**\n- **model**\n\nIf you do not need batch correction, set the method to 'none'.\n\n---\n\n"
    md_content <- paste0(md_content, warning_message)
  }

  # Add DESeq results summary if provided
  if (!is.null(deseq_results)) {
    # Extract total genes with non-zero read count
    total_genes <- sum(!is.na(deseq_results$baseMean) & deseq_results$baseMean > 0)

    # Get thresholds from deseq_results metadata
    metadata_res <- metadata(deseq_results)
    lfc_threshold <- ifelse(!is.null(metadata_res$lfcThreshold), metadata_res$lfcThreshold, 0)
    p_adj_threshold <- ifelse(!is.null(metadata_res$alpha), metadata_res$alpha, 0.1)

    # Extract counts
    res_nonzero <- deseq_results[!is.na(deseq_results$baseMean) & deseq_results$baseMean > 0,]

    # Upregulated genes
    lfc_up <- sum(
      res_nonzero$log2FoldChange > lfc_threshold &
        res_nonzero$padj < p_adj_threshold,
      na.rm = TRUE
    )

    # Downregulated genes
    lfc_down <- sum(
      res_nonzero$log2FoldChange < -lfc_threshold &
        res_nonzero$padj < p_adj_threshold,
      na.rm = TRUE
    )

    # Outliers
    outliers <- sum(
      is.na(res_nonzero$pvalue) & !is.na(res_nonzero$baseMean),
      na.rm = TRUE
    )

    # Low counts (independent filtering)
    low_counts <- sum(
      is.na(res_nonzero$padj) & is.na(res_nonzero$pvalue),
      na.rm = TRUE
    )

    # Mean count threshold (from independent filtering)
    mean_count <- if (!is.null(metadata_res$filterThreshold)) {
      round(metadata_res$filterThreshold, 2)
    } else {
      "-"
    }

    # Calculate percentages
    percent_total <- 100
    percent_up <- if (total_genes > 0) round((lfc_up / total_genes) * 100, 2) else NA
    percent_down <- if (total_genes > 0) round((lfc_down / total_genes) * 100, 2) else NA
    percent_outliers <- if (total_genes > 0) round((outliers / total_genes) * 100, 2) else NA
    percent_low_counts <- if (total_genes > 0) round((low_counts / total_genes) * 100, 2) else NA

    # Create data frame
    summary_df <- data.frame(
      Metric = c(
        "Total Expressed (Non-Zero Read Count)",
        paste0("Upregulated ", "\u2191", " (LFC > ", lfc_threshold, ")"),
        paste0("Downregulated ", "\u2193", " (LFC < ", -lfc_threshold, ")"),
        "Outliers<sup>1</sup>",
        if (mean_count != "-") {
          paste0("Low Counts<sup>2</sup> (Mean Count < ", mean_count, ")")
        } else {
          "Low Counts<sup>2</sup>"
        }
      ),
      `Gene Count` = c(
        total_genes,
        lfc_up,
        lfc_down,
        outliers,
        low_counts
      ),
      `% of Total` = c(
        percent_total,
        percent_up,
        percent_down,
        percent_outliers,
        percent_low_counts
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    # Convert the data frame to markdown format
    # Check if kableExtra is available, otherwise use basic formatting
    if (requireNamespace("kableExtra", quietly = TRUE)) {
      summary_output_md <- paste0(
        "# DESeq2 Summary\n\n",
        kableExtra::kable(
          summary_df,
          format = "html",
          escape = FALSE,
          align = c("l", "c", "c")
        ) %>%
          kableExtra::kable_styling(
            bootstrap_options = c("striped", "hover", "condensed"),
            full_width = FALSE
          ),
        "\n\nArguments of ?DESeq2::results():\n\n",
        "<sup>1</sup> Outliers are genes with high Cook's distance (see [cooksCutoff](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#outlier)).\n\n",
        "<sup>2</sup> Low counts are genes filtered out due to low mean counts (see [independentFiltering](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#independent-filtering-of-results)).\n\n---\n\n"
      )
    } else {
      # Basic markdown table without kableExtra
      table_header <- "| Metric | Gene Count | % of Total |\n|--------|------------|------------|\n"
      table_rows <- apply(summary_df, 1, function(row) {
        paste0("| ", row[1], " | ", row[2], " | ", row[3], " |")
      })
      
      summary_output_md <- paste0(
        "# DESeq2 Summary\n\n",
        table_header,
        paste(table_rows, collapse = "\n"),
        "\n\nArguments of ?DESeq2::results():\n\n",
        "<sup>1</sup> Outliers are genes with high Cook's distance (see cooksCutoff).\n\n",
        "<sup>2</sup> Low counts are genes filtered out due to low mean counts (see independentFiltering).\n\n---\n\n"
      )
    }

    md_content <- paste0(md_content, summary_output_md)
  }

  # Write the content to the output file
  writeLines(md_content, con = output_file)
} 