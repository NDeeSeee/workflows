#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB should be good by default

suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(argparse))
suppressMessages(library(data.table))


set_threads <- function (threads) {
    invisible(capture.output(plan("multiprocess", workers=threads)))
    invisible(capture.output(plan()))
    invisible(capture.output(setDTthreads(threads)))
}


export_rds <- function(data, location){
    tryCatch(
        expr = {
            saveRDS(data, location)
            print(paste("Export data as RDS to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data as RDS to", location, sep=" "))
        }
    )
}


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


load_ctype_data <- function (location) {
    ctype_data_data <- read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    return (ctype_data_data)
}



get_args <- function(){
    parser <- ArgumentParser(description='Assigns cell types to clusters')
    parser$add_argument(
        "--rds",
        help=paste(
            "Path to the RDS file to load Seurat object from. RDS file produced by run_seurat.R script"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--ctype",
        help=paste(
            "Path to the cell types metadata TSV/CSV file with cluster and type columns"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--source",
        help=paste(
            "Column name to select clusters for cell type assignment"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--target",
        help=paste(
            "Column name to store assigned cell types"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./seurat_ctype",
        type="character", default="./seurat_ctype"
    )
    parser$add_argument(
        "--threads",
        help="Threads. Default: 1",
        type="integer", default=1
    )
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}


args <- get_args()
set_threads(args$threads)

print(paste("Loading Seurat data from", args$rds))
seurat_data <- readRDS(args$rds)

print(paste("Loading cell type data from", args$ctype))
ctype_data <- load_ctype_data(args$ctype)

print(paste("Searching for clusters in", args$source, "column. Adding", args$target, "column for cell types"))

seurat_data[[args$target]] <- ctype_data$type[match(
    as.vector(as.character(seurat_data@meta.data[, args$source])),
    as.vector(as.character(ctype_data$cluster)))
]

print("Exporting results")
export_rds(seurat_data, paste(args$output, "_clst_data.rds", sep=""))