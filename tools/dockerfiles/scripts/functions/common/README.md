# Common DESeq2 LRT Analysis Functions

This directory contains common functions used across both Step 1 and Step 2 of the DESeq2 LRT analysis workflow.

## Overview

The DESeq2 LRT (Likelihood Ratio Test) analysis workflow is split into two steps:
1. **Step 1**: Processes raw counts, creates DESeq2 object, and generates contrasts
2. **Step 2**: Analyzes contrasts, generates visualizations, and exports results

To ensure consistency and reduce code duplication, common functionality has been consolidated in this directory.

## Files

- **utilities.R**: Core utility functions including:
  - Metadata validation
  - Batch correction
  - Error handling
  - Helper functions
  - File type handling
  - Step consistency verification

- **visualization.R**: Shared visualization functions:
  - MDS plot creation
  - Common plot styling

- **clustering.R**: Optimized clustering functions:
  - Memory-efficient HOPACH clustering
  - Subset-based clustering for large datasets
  - Scaling methods (z-score, min-max)
  - Various distance metrics (correlation, euclidean)

- **export_functions.R**: Data export functions:
  - GCT file export (for all and filtered data)
  - DESeq2 result export
  - Contrast information export
  - Analysis parameter export
  - CLS export for GSEA
  - Visualization export

## Usage

These functions are sourced at the beginning of both the Step 1 and Step 2 scripts:

```r
# Source common functions
source("functions/common/utilities.R")
source("functions/common/visualization.R")
source("functions/common/clustering.R")
source("functions/common/export_functions.R")
```

## Design Principles

1. **Memory Efficiency**: Functions are optimized for large datasets with memory management
2. **Error Handling**: Comprehensive error checking with informative messages
3. **Consistency**: Ensures consistent behavior between Step 1 and Step 2
4. **Flexibility**: Functions accept various parameters to adapt to different use cases

## Consolidation

These common functions were consolidated from previously separate implementations in the Step 1 and Step 2 directories to:
- Reduce code duplication
- Ensure consistent behavior
- Simplify maintenance
- Improve organization of the codebase 