# install_packages.R

# List of CRAN packages to install
cran_packages <- c(
  "locfit", "tidyverse", "patchwork", "argparse", "data.table",
  "pheatmap", "ggrepel", "htmlwidgets", "devtools", "kableExtra"
)

# List of Bioconductor packages to install
bioc_packages <- c(
  "BiocParallel", "DESeq2", "limma", "EnhancedVolcano", "hopach", "cmapR"
)

# Function to install and check CRAN packages
install_and_check_cran <- function(packages) {
  install.packages(packages, repos = "https://cloud.r-project.org")
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("Package", pkg, "failed to install."))
    }
  }
}

# Function to install and check Bioconductor packages
install_and_check_bioc <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(version = "3.19", update = FALSE, ask = FALSE)
  BiocManager::install(packages, update = TRUE, ask = FALSE)
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("Package", pkg, "failed to install."))
    }
  }
}

# Install CRAN packages
install_and_check_cran(cran_packages)

# Install Bioconductor packages
install_and_check_bioc(bioc_packages)

# Install additional packages from GitHub
if (!require("devtools", character.only = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
devtools::install_github("hasaru-k/GlimmaV2")
if (!require("GlimmaV2", character.only = TRUE)) {
  stop("Package GlimmaV2 failed to install from GitHub.")
}

cat("All packages installed successfully.\n")