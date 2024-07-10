#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(argparse))
suppressMessages(library(DropletUtils))


get_args <- function(){
    parser <- ArgumentParser(description="Splits feature-barcode matrix produced by Cell Ranger ARC into RNA and ATAC matrices")
    parser$add_argument(
        "--mex",
        help=paste(
            "Path to the folder with feature-barcode matrix from Cell Ranger ARC Count/Aggregate",
            "experiment in MEX format. The rows consist of all the genes and peaks concatenated",
            "together and the columns are restricted to those barcodes that are identified as cells."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./sc",
        type="character", default="./sc"
    )
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    print(args)
    return (args)
}

args <- get_args()

print(paste("Loading counts from", args$mex))
all_counts <- Seurat::Read10X(data.dir=args$mex)
rna_counts <- all_counts[["Gene Expression"]]
atac_counts <- all_counts[["Peaks"]]

rna_outputs_location <- paste(args$output, "rna", sep="_")
print(paste("Exporting RNA counts to", rna_outputs_location))
DropletUtils:::write10xCounts(
    path=rna_outputs_location,
    x=rna_counts,
    type="sparse",
    version="3",
    overwrite=TRUE
)

atac_outputs_location <- paste(args$output, "atac", sep="_")
print(paste("Exporting ATAC counts to", atac_outputs_location))
DropletUtils:::write10xCounts(
    path=atac_outputs_location,
    x=atac_counts,
    type="sparse",
    version="3",
    overwrite=TRUE
)