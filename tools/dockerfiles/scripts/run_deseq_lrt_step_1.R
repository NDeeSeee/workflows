#!/usr/bin/env Rscript
options(warn = -1)
options("width" = 400)

suppressMessages(library(argparse)) # For argument parsing
suppressMessages(library(BiocParallel)) # For parallel processing
suppressMessages(library(pheatmap)) # For heatmap
suppressMessages(library(DESeq2)) # For DESeq2
suppressMessages(library(tidyverse)) # For data manipulation
suppressMessages(library(Glimma)) # For glimmaMDS
suppressMessages(library(cmapR)) # For write_gct
suppressMessages(library(sva)) # For ComBat_seq
suppressMessages(library(limma)) # For removeBatchEffect
suppressMessages(library(hopach)) # For clustering
suppressMessages(library(rlang)) # For traceback handling

options(rlang_backtrace_on_error = "full")

mutate <- dplyr::mutate
filter <- dplyr::filter
group_by <- dplyr::group_by
slice <- dplyr::slice
rename <- dplyr::rename
select <- dplyr::select
`%>%` <- magrittr::`%>%`
`%in%` <- base::`%in%`

###### v0.0.5 ######
#
# Changes:
# - 'combatseq' batch correction is applied in Step 1 if specified.
# - 'limmaremovebatcheffect' is noted and passed to Step 2 for application.
# - Batch correction options are handled appropriately.
# - Necessary information is saved for Step 2.
#
###### v0.0.4 ######
#
# Changes:
# - Adjusted batch correction to be optional.
# - Batch correction with limma is only applied if specified.
#
###### v0.0.3 ######
#
# Changes:
# - Reads preprocessed data from RDS files.
# - Accepts `contrast_indices`, `fdr`, `lfcThreshold`, and `regulation` as inputs.
# - Applies batch correction using `limma::removeBatchEffect` if specified.
# - Uses `results` function with appropriate parameters.
#
###### v0.0.2 ######
#
# Changes:
# - Save DESeq2 dataset objects (dds) and LRT results (dsq_lrt_res) into RDS files.
# - Save expression data (expression_data_df) into an RDS file.
# - Save metadata (metadata_df) into an RDS file.
# - Adjust outputs to include these RDS files.
#
###### v0.0.1 ######
#
# Note: at least two biological replicates are required for every compared category.
#
# All input CSV/TSV files should have the following header (case-sensitive)
# <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
# <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV
#
# Format of the input files is identified based on file's extension
# *.csv - CSV
# *.tsv - TSV
# Otherwise used CSV by default
#
# The output file's rows order corresponds to the rows order of the first CSV/TSV file.
# Output file is always saved in TSV format
#
# Output file includes only intersected rows from all input files. Intersected by
# RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand
#
# Additionally we calculate -LOG10(pval) and -LOG10(padj)
#
# Example of CSV metadata file set with --meta
#
# ,time,condition
# DH1,day5,WT
# DH2,day5,KO
# DH3,day7,WT
# DH4,day7,KO
# DH5,day7,KO
#
# where time, condition, day5, day7, WT, KO should be a single words (without spaces)
# and DH1, DH2, DH3, DH4, DH5 correspond to the --names (spaces are allowed)
#
# --contrast should be set based on your metadata file in a form of Factor Numerator Denominator
# where Factor      - columns name from metadata file
#       Numerator   - category from metadata file to be used as numerator in fold change calculation
#       Denominator - category from metadata file to be used as denominator in fold change calculation
# for example condition WT KO
# if --contrast is set as a single string "condition WT KO" then is will be splitted by space
#

###### run_deseq_lrt_step_1.R ######

READ_COL <- "TotalReads"
RPKM_COL <- "Rpkm"
INTERSECT_BY <- c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand")

get_file_type <- function(filename) {
  ext <- tools::file_ext(filename)
  separator <- ","
  if (ext == "tsv") {
    separator <- "\t"
  }
  return(separator)
}

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
      collected_isoforms <- merge(collected_isoforms,
        isoforms,
        by = intersect_by,
        sort = FALSE
      )
    }
  }
  print(paste0("Number of rows common for all loaded files ", nrow(collected_isoforms)))

  return(collected_isoforms)
}

assert_args <- function(args) {
  if (length(args$input) != length(args$name)) {
    print("Exiting: --input and --name have different number of values")
    quit(
      save = "no",
      status = 1,
      runLast = FALSE
    )
  }
  tryCatch(
    expr = {
      # Try to load design formula
      design_formula <- as.formula(args$design)
    },
    error = function(e) {
      print(paste0("Exiting: failed to load --design ", args$design, " as formula"))
      quit(
        save = "no",
        status = 1,
        runLast = FALSE
      )
    }
  )
  tryCatch(
    expr = {
      # Try to load reduced formula
      reduced_formula <- as.formula(args$reduced)
    },
    error = function(e) {
      print(paste0("Exiting: failed to load --reduced ", args$reduced, " as formula"))
      quit(
        save = "no",
        status = 1,
        runLast = FALSE
      )
    }
  )
  return(args)
}

# Function to apply batch correction using ComBat_seq
apply_combat_seq <- function(count_data, batch_info, design, column_data) {
  print("Applying ComBat_seq batch correction")
  corrected_counts <- sva::ComBat_seq(
    counts = as.matrix(count_data),
    batch = batch_info,
    group = NULL,
    covar_mod = model.matrix(design, data = column_data)
  )
  return(corrected_counts)
}

get_args <- function() {
  parser <- ArgumentParser(description = "Run DESeq2 for multi-factor analysis using LRT (likelihood ratio or chi-squared test)")
  parser$add_argument(
    "-i",
    "--input",
    help = "Grouped by gene / TSS/ isoform expression files, formatted as CSV/TSV",
    type = "character",
    required = "True",
    nargs = "+"
  )
  parser$add_argument(
    "-n",
    "--name",
    help = "Unique names for input files, no special characters, spaces are allowed. Number and order corresponds to --input",
    type = "character",
    required = "True",
    nargs = "+"
  )
  parser$add_argument(
    "-m",
    "--meta",
    help = "Metadata file to describe relation between samples, where first column corresponds to --name, formatted as CSV/TSV",
    type = "character",
    required = "True"
  )
  parser$add_argument(
    "-d",
    "--design",
    help = "Design formula. Should start with ~, like ~condition+celltype+condition:celltype. See DESeq2 manual for details",
    type = "character",
    required = "True"
  )
  parser$add_argument(
    "-r",
    "--reduced",
    help = "Reduced formula to compare against with the term(s) of interest removed. Should start with ~. See DESeq2 manual for details",
    type = "character",
    required = "True"
  )
  parser$add_argument(
    "--batchcorrection",
    help = paste(
      "Specifies the batch correction method to be applied.",
      "- 'combatseq' applies ComBat_seq at the beginning of the analysis, removing batch effects from the counts before differential expression analysis.",
      "- 'limmaremovebatcheffect' applies removeBatchEffect from the limma package after differential expression analysis.",
      "- Default: none"
    ),
    type = "character",
    choices = c("none", "combatseq", "limmaremovebatcheffect"),
    default = "none"
  )
  parser$add_argument(
    "--fdr",
    help = paste(
      "In the exploratory visualization part of the analysis output only features",
      "with adjusted p-value (FDR) not bigger than this value. Also the significance",
      "cutoff used for optimizing the independent filtering. Default: 0.1."
    ),
    type = "double",
    default = 0.1
  )
  parser$add_argument(
    "--lfcthreshold",
    help = paste(
      "Log2 fold change threshold for determining significant differential expression.",
      "Genes with absolute log2 fold change greater than this threshold will be considered.",
      "Default: 0.59 (about 1.5 fold change)"
    ),
    type = "double",
    default = 0.59
  )
  parser$add_argument(
    "--use_lfc_thresh",
    help = paste(
      "Flag to indicate whether to use lfcthreshold as the null hypothesis value in the results function call.",
      "If TRUE, lfcthreshold is used in the hypothesis test (i.e., genes are tested against this threshold).",
      "If FALSE, the null hypothesis is set to 0, and lfcthreshold is used only as a downstream filter.",
      "Default: FALSE"
    ),
    action = "store_true",
    default = FALSE
  )
  parser$add_argument(
    "--cluster",
    help    = paste(
      "Hopach clustering method to be run on normalized read counts for the",
      "exploratory visualization part of the analysis. Default: do not run",
      "clustering"
    ),
    type    = "character",
    choices = c("row", "column", "both", "none"),
    default = "none"
  )
  parser$add_argument(
    "-o",
    "--output",
    help = "Output prefix for generated files",
    type = "character",
    default = "./deseq"
  )
  parser$add_argument(
    "-p",
    "--threads",
    help = "Threads number",
    type = "integer",
    default = 1
  )
  parser$add_argument(
    "--test_mode",
    help = "Run for test, only first 500 rows",
    action = "store_true",
    default = FALSE
  )
  args <- assert_args(parser$parse_args(commandArgs(trailingOnly = TRUE)))
  return(args)
}

generate_lrt_md <- function(deseq_results, full_formula, reduced_formula, output_file, alpha = args$fdr,
                            batch_warning = NULL) {
  # Initialize the markdown content
  md_content <- ""

  # Add batch warning if present
  if (!is.null(batch_warning)) {
    md_content <- paste0(md_content, "# **Warning**\n\n", batch_warning, "\n\n---\n\n")
  }

  # Add DESeq LRT results summary based on padj value only
  if (!is.null(deseq_results)) {
    # Start summarizing the LRT results
    md_content <- paste0(md_content, "# Likelihood Ratio Test (LRT) Results\n\n---\n\n")

    # Describe the full and reduced formulas
    md_content <- paste0(
      md_content,
      "Based on your **full formula**: `", full_formula, "` and **reduced formula**: `", reduced_formula, "`, ",
      "this LRT analysis tests whether removing the interaction term (or terms) significantly affects gene expression. ",
      "The test uses only the **FDR adjusted p-value** (padj) to determine significance, as Log Fold Change (LFC) is irrelevant in the context of LRT.\n\n"
    )

    # Calculate the number of significant genes based on padj value
    significant_genes <- sum(deseq_results$padj < alpha, na.rm = TRUE)
    total_genes <- sum(!is.na(deseq_results$padj))

    # Extract the outliers and low counts from the DESeq results summary
    summary_output <- capture.output(summary(deseq_results))

    outliers <- gsub(".*: ", "", summary_output[6]) # Outliers
    low_counts <- gsub(".*: ", "", summary_output[7]) # Low counts (independent filtering)
    mean_count <- gsub("[^0-9]", "", summary_output[8])

    # Add a summary of significant genes
    lrt_summary <- paste0(
      "### Results Summary\n\n",
      "From this LRT analysis, **", significant_genes, " genes** (out of ", total_genes, " tested) are identified as significant with a padj value < ", alpha, ".\n\n"
    )

    # Add information about outliers and low counts
    lrt_summary <- paste0(
      lrt_summary,
      "**Outliers**<sup>1</sup>: ", outliers, " of genes were detected as outliers and excluded from analysis.\n\n",
      "**Low counts**<sup>2</sup>: ", low_counts, " of genes were removed due to low counts (mean <", mean_count, ") and independent filtering.\n\n",
      "Arguments of ?DESeq2::results():   \n<sup>1</sup> - see 'cooksCutoff',\n<sup>2</sup> - see 'independentFiltering'\n\n"
    )

    # Add this summary to the markdown content
    md_content <- paste0(md_content, lrt_summary)

    # Add explanation for the next steps
    next_steps <- paste0(
      "---\n\n",
      "### Next Steps\n\n",
      "If the number of significant genes is substantial, consider including the interaction term in your design formula ",
      "for a more detailed analysis of differential gene expression.\n\n",
      "For further insights and to explore detailed contrasts using the Wald test for the complex design formula, ",
      "please visit the **Complex Interaction Analysis** tab for more information.\n\n"
    )

    md_content <- paste0(md_content, next_steps)
  }

  # Write the content to the output file
  writeLines(md_content, con = output_file)
}

# Function to generate main effect contrasts with different reference levels
generate_main_effect_contrasts <- function(dds, factors, factor_levels) {
  contrasts <- list()

  for (factor in factors) {
    other_factors <- setdiff(factors, factor)

    for (other_factor in other_factors) {
      other_levels <- factor_levels[[other_factor]]

      for (other_level in other_levels) {
        # dds_subset <- dds[colData(dds)[[other_factor]] == other_level, ]
        dds_subset <- dds

        for (ref_level in factor_levels[[factor]]) {
          colData(dds_subset)[[factor]] <- relevel(colData(dds_subset)[[factor]], ref = ref_level)
          colData(dds_subset)[[other_factor]] <- relevel(colData(dds_subset)[[other_factor]], ref = other_level)

          levels <- factor_levels[[factor]]

          for (level in levels) {
            if (level != ref_level) {
              specificity_group <- paste(other_factor, other_level, sep = "_")
              contrast <- paste(factor, level, "vs", ref_level, sep = "_")

              dds_subset <- DESeq(dds_subset, test = "Wald")

              contrast_res <- results(dds_subset,
                name = contrast,
                alpha = args$fdr,
                lfcThreshold = ifelse(args$use_lfc_thresh, args$lfcthreshold, 0),
                independentFiltering = TRUE
              )

              significant_genes <- sum(contrast_res$padj < args$fdr, na.rm = TRUE)

              contrasts <- append(contrasts, list(list(
                effect_type = "main",
                specificity_group = specificity_group,
                numerator = level,
                denominator = ref_level,
                contrast = contrast,
                subset = dds_subset,
                contrast_res = contrast_res,
                significant_genes = significant_genes
              )))
            }
          }
        }
      }
    }
  }

  return(contrasts)
}

# Function to dynamically extract the factor and level names from interaction terms
extract_factors_and_levels <- function(term) {
  parts <- strsplit(term, "\\.")[[1]]
  factor1 <- sub("[0-9A-Z_]+$", "", parts[1])
  factor2 <- sub("[0-9A-Z_]+$", "", parts[2])
  level1  <- sub("^.*?([0-9A-Z_]+$)", "\\1", parts[1])
  level2  <- sub("^.*?([0-9A-Z_]+$)", "\\1", parts[2])
  list(factor1 = factor1, level1 = level1, factor2 = factor2, level2 = level2)
}

# Function to generate interaction effect contrasts
generate_interaction_effect_contrasts <- function(dds) {
  contrasts <- list()
  interaction_names <- resultsNames(dds)

  # Identify interaction terms in the results names
  interaction_terms <- grep("\\.", interaction_names, value = TRUE)

  for (interaction in interaction_terms) {
    factors_levels <- extract_factors_and_levels(interaction)
    factor1 <- factors_levels$factor1
    level1 <- factors_levels$level1
    factor2 <- factors_levels$factor2
    level2 <- factors_levels$level2

    levels1 <- levels(colData(dds)[[factor1]])
    levels2 <- levels(colData(dds)[[factor2]])

    # Generate contrasts for factor1 vs factor2
    for (ref_level1 in levels1) {
      for (ref_level2 in levels2) {
        if ((ref_level1 != level1 || ref_level2 != level2) && (ref_level2 != level2)) {
          specificity_group <- paste(factor2, level2, "vs", ref_level2, sep = "_")
          numerator <- paste0(factor1, level1)
          denominator <- paste0(factor1, ref_level1)

          if (numerator != denominator && specificity_group != paste(factor2, ref_level2, "vs", ref_level2, sep = "_")) {
            # dds_subset <- dds[colData(dds)[[factor2]] == ref_level2, ]
            dds_subset <- dds
            colData(dds_subset)[[factor1]] <- relevel(colData(dds_subset)[[factor1]], ref = ref_level1)

            dds_subset <- DESeq(dds_subset, test = "Wald")

            contrast_res <- results(dds_subset,
              name = interaction,
              alpha = args$fdr,
              lfcThreshold = ifelse(args$use_lfc_thresh, args$lfcthreshold, 0),
              independentFiltering = TRUE
            )

            significant_genes <- sum(contrast_res$padj < args$fdr, na.rm = TRUE)

            contrasts <- append(contrasts, list(list(
              effect_type = "interaction",
              specificity_group = specificity_group,
              numerator = numerator,
              denominator = denominator,
              contrast = interaction,
              subset = dds_subset,
              contrast_res = contrast_res,
              significant_genes = significant_genes
            )))
          }
        }
      }
    }

    # Generate contrasts for factor2 vs factor1
    for (ref_level2 in levels2) {
      for (ref_level1 in levels1) {
        if ((ref_level2 != level2 || ref_level1 != level1) && (ref_level1 != level1)) {
          specificity_group <- paste(factor1, level1, "vs", ref_level1, sep = "_")
          numerator <- paste0(factor2, level2)
          denominator <- paste0(factor2, ref_level2)

          if (numerator != denominator && specificity_group != paste(factor1, ref_level1, "vs", ref_level1, sep = "_")) {
            # dds_subset <- dds[colData(dds)[[factor1]] == ref_level1, ]
            dds_subset <- dds
            colData(dds_subset)[[factor2]] <- relevel(colData(dds_subset)[[factor2]], ref = ref_level2)

            dds_subset <- DESeq(dds_subset, test = "Wald")

            contrast_res <- results(dds_subset,
              name = interaction,
              alpha = args$fdr,
              lfcThreshold = ifelse(args$use_lfc_thresh, args$lfcthreshold, 0),
              independentFiltering = TRUE
            )

            significant_genes <- sum(contrast_res$padj < args$fdr, na.rm = TRUE)

            contrasts <- append(contrasts, list(list(
              effect_type = "interaction",
              specificity_group = specificity_group,
              numerator = numerator,
              denominator = denominator,
              contrast = interaction,
              subset = dds_subset,
              contrast_res = contrast_res,
              significant_genes = significant_genes
            )))
          }
        }
      }
    }
  }

  return(contrasts)
}

# Main function to generate the dataframe with all possible contrasts
generate_contrasts <- function(dds) {
  # Extract the design formula and model matrix
  design_formula <- design(dds)

  # Get the levels of each factor in the design
  factors <- all.vars(design_formula)
  factors <- factors[!str_detect(factors, "batch")]
  print("Analyzing the following factors:")
  print(factors)

  factor_levels <- lapply(factors, function(f) levels(colData(dds)[[f]]))
  names(factor_levels) <- factors
  print("Analyzing the following factor levels:")
  print(factor_levels)

  # Generate main and interaction effect contrasts
  main_contrasts <- generate_main_effect_contrasts(dds, factors, factor_levels)
  interaction_contrasts <- generate_interaction_effect_contrasts(dds)
  all_contrasts <- c(main_contrasts, interaction_contrasts)

  all_contrasts <- purrr::list_modify(all_contrasts, expression_data_df = expression_data_df)

  print(paste("Exporting contrasts list to", paste0(args$output, "_contrasts.rds"), sep = " "))
  saveRDS(all_contrasts, file = paste0(args$output, "_contrasts.rds"))

  # Create a dataframe to store the contrasts and number of significant genes
  contrast_df <- data.frame(
    effect = character(),
    specificity_group = character(),
    contrast = character(),
    numerator = character(),
    denominator = character(),
    significant_genes = integer(),
    stringsAsFactors = FALSE
  )

  # Run DESeq2 for each contrast and store the results
  count <- 1
  for (contrast in all_contrasts) {
    print(count)

    contrast_df <- rbind(contrast_df, data.frame(
      effect = contrast$effect_type,
      contrast = contrast$contrast,
      specificity_group = contrast$specificity_group,
      numerator = contrast$numerator,
      denominator = contrast$denominator,
      significant_genes = contrast$significant_genes,
      stringsAsFactors = FALSE
    ))
    count <- count + 1
  }

  # Remove duplicate contrasts
  contrast_df <- contrast_df %>%
    group_by(specificity_group) %>%
    distinct(contrast, .keep_all = TRUE) %>%
    arrange(specificity_group) %>%
    ungroup() %>%
    ungroup() %>%
    mutate(contrast_number = row_number()) %>%
    select(contrast_number, everything()) %>%
    arrange(desc(significant_genes))

  # Sort the dataframe by the number of significant genes
  return(contrast_df)
}

# Function to clean sample names
clean_sample_names <- function(names) {
  names <- trimws(names) # Remove leading/trailing whitespace
  names <- tolower(names) # Convert to lowercase
  names <- gsub("\\s+", "_", names) # Replace spaces with underscores
  names <- gsub("[^[:alnum:]_]", "", names) # Remove non-alphanumeric characters except underscores
  names <- gsub("^_+|_+$", "", names) # Remove leading/trailing underscores
  return(names)
}

# Function to export MDS plot
export_mds_html_plot <- function(norm_counts_data, location) {
  tryCatch(
    expr  = {
      htmlwidgets::saveWidget(
        Glimma::glimmaMDS(
          x      = norm_counts_data,
          groups = as.data.frame(metadata_df),
          labels = rownames(metadata_df)
        ),
        file = location
      )
    },
    error = function(e) {
      print(paste0("Failed to export MDS plot to ", location, " with error - ", e))
    }
  )
}

# Function to export normalized counts and filtered counts to GCT format
export_gct_data <- function(normCounts, row_metadata, col_metadata, output_prefix) {
  tryCatch(
    expr = {
      # Use the provided row_metadata directly
      # Prepare column metadata
      col_metadata <- col_metadata %>% mutate_all(as.vector)

      # Create GCT object
      gct_data <- new("GCT",
                      mat   = normCounts,
                      rdesc = row_metadata,
                      cdesc = col_metadata)

      # Write GCT files
      cmapR::write_gct(
        ds    = gct_data,
        ofile = paste0(output_prefix, "_counts_all.gct"),
        appenddim = FALSE
      )
      print(paste("Exporting GCT data to", paste0(output_prefix, "_counts_all.gct")))

      # Filter rows by padj <= FDR if padj is available
      if ('padj' %in% colnames(row_metadata)) {
        row_metadata_filtered <- row_metadata %>% dplyr::filter(padj <= args$fdr)
      } else {
        row_metadata_filtered <- row_metadata
      }

      # Check if any genes remain after filtering
      if (nrow(row_metadata_filtered) == 0) {
        warning(paste("No genes passed the FDR threshold of", args$fdr, "for", output_prefix))
      }

      # Subset normCounts to include only the filtered rows
      filtered_normCounts <- normCounts[rownames(row_metadata_filtered),]

      # Create filtered GCT object
      gct_data_filtered <- new("GCT",
                               mat   = filtered_normCounts,
                               rdesc = row_metadata_filtered,
                               cdesc = col_metadata)

      # Write filtered GCT file
      cmapR::write_gct(
        ds    = gct_data_filtered,
        ofile = paste0(output_prefix, "_counts_filtered.gct"),
        appenddim = FALSE
      )
      print(paste("Exporting GCT data to", paste0(output_prefix, "_counts_filtered.gct")))
    },
    error = function(e) {
      print(paste("Failed to export GCT data to", output_prefix, "with error -", e$message))
    }
  )
}

export_charts <- function(res, annotated_expression_df, column_data, normCounts, output, args) {
  # Export MDS plot
  print("Exporting MDS plot")
  export_mds_html_plot(normCounts, paste0(output, "_mds_plot.html"))

  print("Exporting MDS plot with batch correction")
  if (!is.null(corrected_counts)) {
    export_mds_html_plot(corrected_counts, paste0(output, "_mds_plot_corrected.html"))
  }

  clustered_data <- cluster_and_reorder(normCounts, column_data, annotated_expression_df, args)
  # Export GCT data
  export_gct_data(clustered_data$normCounts, clustered_data$row_metadata, clustered_data$col_metadata, output)
}

# Function to generate clusters
get_clustered_data <- function(expression_data, transpose = FALSE) {
  if (transpose) {
    print("Transposing expression data")
    expression_data <- t(expression_data)
  }

  # Apply scaling per row and maintain the original orientation
  expression_data <- t(apply(
    expression_data,
    1,
    FUN = function(x) {
      scale_min_max(x)
    }
  ))

  if (transpose) {
    print("Transposing expression data back")
    expression_data <- t(expression_data)
  }

  print("Running HOPACH")
  hopach_results <- hopach::hopach(expression_data, verbose = TRUE)

  print("Parsing cluster labels")
  clusters <- data.frame(label = hopach_results$clustering$labels)
  clusters <- clusters %>% mutate(HCL = paste0("c", label))

  return(list(
    order    = as.vector(hopach_results$clustering$order),
    expression = expression_data,
    clusters = clusters
  ))
}

# Function for min-max scaling
scale_min_max <- function(x,
                          min_range = -2,
                          max_range = 2) {
  min_val  <- min(x)
  max_val  <- max(x)
  scaled_x <-
    (x - min_val) / (max_val - min_val) * (max_range - min_range) + min_range
  return(scaled_x)
}

# Function to cluster data and re-order based on clustering results
cluster_and_reorder <- function(normCounts, col_metadata, row_metadata, args) {
  if (args$cluster != "none") { # Only cluster if not "none"
    if (args$cluster == "column" || args$cluster == "both") {
      print("Clustering filtered read counts by columns")
      clustered_data <- get_clustered_data(expression_data = normCounts, transpose = TRUE)
      col_metadata <- cbind(col_metadata, clustered_data$clusters)
      col_metadata <- col_metadata[clustered_data$order,]
      normCounts   <- normCounts[, clustered_data$order]
      print("Reordered samples")
      print(col_metadata)
    }
    if (args$cluster == "row" || args$cluster == "both") {
      print("Clustering filtered normalized read counts by rows")
      clustered_data <- get_clustered_data(expression_data = normCounts, transpose = FALSE)
      normCounts   <- clustered_data$expression
      row_metadata <- cbind(row_metadata, clustered_data$clusters)
      row_metadata <- row_metadata[clustered_data$order,]
      normCounts   <- normCounts[clustered_data$order,]
      print("Reordered features")
      print(head(row_metadata))
    }
  } else {
    print("Clustering skipped as per 'none' option.")
  }
  return(list(normCounts = normCounts, col_metadata = col_metadata, row_metadata = row_metadata))
}

# Parse arguments
args <- get_args()

# Set threads
register(MulticoreParam(args$threads))

# Load metadata
metadata_df <- read.table(
  args$meta,
  sep = get_file_type(args$meta),
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = 1
)
print(paste("Load metadata from", args$meta, sep = " "))
print(metadata_df)

# Load design formula
design_formula <- as.formula(args$design)
print("Load design formula")
print(design_formula)

# Load reduced formula
reduced_formula <- as.formula(args$reduced)
print("Load reduced formula")
print(reduced_formula)

# Clean sample names
args$name <- clean_sample_names(args$name)
colnames(metadata_df) <- clean_sample_names(colnames(metadata_df))
rownames(metadata_df) <- clean_sample_names(rownames(metadata_df))

print("Processed args$name:")
print(args$name)

print("Processed rownames(metadata_df):")
print(rownames(metadata_df))

# Load expression data
expression_data_df <- load_expression_data(args$input, args$name, READ_COL, RPKM_COL, INTERSECT_BY)
print("Expression data")
print(head(expression_data_df))

# Select all columns with read counts data, reorder them based on the row names from metadata_df
read_counts_columns <- grep(
  paste(READ_COL, sep = ""),
  colnames(expression_data_df),
  value = TRUE,
  ignore.case = TRUE
)

# TODO: double-check if we should use GeneId or RefseqId here
read_counts_data_df <- expression_data_df %>%
  dplyr::mutate_at("GeneId", toupper) %>%
  dplyr::distinct(GeneId, .keep_all = TRUE) %>% # to prevent from failing when input files are not grouped by GeneId
  remove_rownames() %>%
  column_to_rownames("GeneId")

read_counts_data_df <- read_counts_data_df[read_counts_columns]

# Clean column names of read_counts_data_df
colnames(read_counts_data_df) <- lapply(colnames(read_counts_data_df), function(s) {
  paste(head(unlist(strsplit(s, " ", fixed = TRUE)), -1), collapse = " ")
})
colnames(read_counts_data_df) <- clean_sample_names(colnames(read_counts_data_df))

print("Processed colnames(read_counts_data_df):")
print(colnames(read_counts_data_df))

print("Read counts data")
print(head(read_counts_data_df))

# Check for mismatches in sample names
samples_in_metadata_not_in_data <- setdiff(rownames(metadata_df), colnames(read_counts_data_df))
samples_in_data_not_in_metadata <- setdiff(colnames(read_counts_data_df), rownames(metadata_df))

if (length(samples_in_metadata_not_in_data) > 0) {
  print("Error: The following samples are in the metadata but not in the expression data:")
  print(samples_in_metadata_not_in_data)
  stop("Sample names mismatch between metadata and expression data.")
}

if (length(samples_in_data_not_in_metadata) > 0) {
  print("Error: The following samples are in the expression data but not in the metadata:")
  print(samples_in_data_not_in_metadata)
  stop("Sample names mismatch between expression data and metadata.")
}

corrected_counts <- NULL

# Reorder read counts data based on metadata
read_counts_data_df <- read_counts_data_df[, rownames(metadata_df)]

# Initialize batch warning message
batch_warning <- NULL

countData <- read_counts_data_df

# Check for batch correction and handle accordingly
if (args$batchcorrection != "none") {
  if ("batch" %in% colnames(metadata_df)) {
    # Ensure 'batch' is a factor
    metadata_df$batch <- as.factor(metadata_df$batch)

    if (args$batchcorrection == "combatseq") {
      # Case 2: Apply ComBat_seq batch correction to raw counts
      print("Applying ComBat_seq batch correction")
      batch_info       <- metadata_df$batch
      corrected_counts <- apply_combat_seq(countData, batch_info, design_formula, metadata_df)
      countData        <- corrected_counts
      # Design formula remains as provided
    } else if (args$batchcorrection == "limmaremovebatcheffect") {
      # Case 3: Include 'batch' in the design formula
      print("Including 'batch' in the design formula for limma batch correction")
      # Remove '~' from the original design formula string
      original_design <- substring(args$design, 2)
      # Create new design formula with 'batch' included
      design_formula  <- as.formula(paste("~ batch +", original_design))
      # Note: We do not modify counts at this stage
    }
  } else {
    batch_warning <- "You provided batch correction but there's no 'batch' column in your metadata."
    print(batch_warning)
    # Proceed without batch correction
    # countData remains as read_counts_data_df
  }
} else {
  # Case 1: No batch correction
  # countData remains as read_counts_data_df
  # Design formula remains as provided
  print("No batch correction provided; proceeding with default settings")
}

# Create DESeq2 dataset
dse <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = metadata_df,
  design = design_formula
)

# Extract the raw counts matrix
cts <- counts(dse)  # a matrix of counts, rows = genes, columns = samples

# Sum counts per sample (sum columns to get total reads per sample)
sample_sums <- colSums(cts)

# Create a data frame for plotting
plot_df <- data.frame(
  Sample = names(sample_sums),
  Count  = sample_sums
)

# Now 'plot_df' has three columns: Sample, Count, and Statistic
# You can plot them using ggplot2, similar to the previous logic

stat_barchart <- ggplot(plot_df, aes(x = Sample, y = Count, fill = "blue")) +
  geom_col(position = "dodge") +
  labs(title = "Total Reads by Sample",
       x     = "Sample",
       y     = "Count",
       fill  = "blue") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("alignment_stats_barchart.png",
       plot   = stat_barchart,
       dpi    = 400,
       height = 3,
       width  = max(10, length(sample_sums) * 0.5),
       units  = "in")

print("Run DESeq2 using Wald")
dsq_wald <- DESeq(
  dse,
  test = "Wald",
  quiet = FALSE,
  parallel = TRUE
)

# Obtain normalized counts for downstream analyses
if (args$batchcorrection == "limmaremovebatcheffect" && "batch" %in% colnames(metadata_df)) {
  # Case 3: After DESeq2, apply rlog transformation and remove batch effects using limma
  print("Applying rlog transformation and limma batch effect removal")
  rlog_transformed <- rlog(dsq_wald, blind = FALSE)
  rlog_counts      <- assay(rlog_transformed)

  # Prepare design matrix without 'batch' for removeBatchEffect
  design_formula <- as.formula(paste0("~", str_remove(as.character(dsq_wald@design)[2], " \\+ batch")))

  design_matrix <- model.matrix(design_formula, data = metadata_df)

  # Apply removeBatchEffect
  corrected_counts <- limma::removeBatchEffect(rlog_counts, batch = metadata_df$batch, design = design_matrix)
  normCounts       <- corrected_counts
} else {
  # Case 1 and Case 2: Use rlog-transformed counts without additional batch correction
  print("Applying rlog transformation without additional batch effect removal")
  rlog_transformed <- rlog(dsq_wald, blind = FALSE)
  normCounts       <- assay(rlog_transformed)
}

print("Run DESeq2 using LRT")
dsq_lrt <- DESeq(
  dse,
  test = "LRT",
  reduced = reduced_formula,
  quiet = FALSE,
  parallel = TRUE
)

print("Generate contrasts")
contrast_df <- generate_contrasts(dsq_wald)

dsq_lrt_res <- results(dsq_lrt,
  alpha = args$fdr,
  independentFiltering = TRUE
)

print("LRT Results description")
print(mcols(dsq_lrt_res))

# Export results to TSV file
annotated_expression_df <- expression_data_df %>%
  bind_cols(as.data.frame(dsq_lrt_res) %>% select(pvalue, padj)) %>%
  mutate(
    `-LOG10(pval)` = -log10(as.numeric(pvalue)),
    `-LOG10(padj)` = -log10(as.numeric(padj))
  )

contrasts_filename <- paste(args$output, "_contrasts_table.tsv", sep = "")
write.table(
  contrast_df,
  file = contrasts_filename,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

lrt_report_filename <- paste0(args$output, "_lrt_result.md")
summary(dsq_lrt_res)
generate_lrt_md(dsq_lrt_res, args$design, args$reduced, lrt_report_filename, alpha = args$fdr, batch_warning = batch_warning)
print(paste("Export LRT markdown report to", lrt_report_filename, sep = " "))

results_filename <- paste0(args$output, "_gene_exp_table.tsv")
write.table(
  annotated_expression_df,
  file = results_filename,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
print(paste("Export results to", results_filename, sep = " "))
graphics.off()

export_charts(dsq_lrt_res, annotated_expression_df, metadata_df, normCounts, args$output, args)