#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB should be good by default

suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(argparse))


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


export_data <- function(data, location, row_names=FALSE, col_names=TRUE, quote=FALSE){
    tryCatch(
        expr = {
            write.table(
                data,
                file=location,
                sep=get_file_type(location),
                row.names=row_names,
                col.names=col_names,
                quote=quote
            )
            print(paste("Export data to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data to", location, sep=" "))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description="Extracts cells, UMAP, and clusters from the Seurat RDS file")
    parser$add_argument(
        "--rds",
        help=paste(
            "Path to the RDS file to load Seurat object from.",
            "RDS file can be produced by run_seurat.R or assign_cell_types.R scripts"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--source",
        help=paste(
            "Column name to extract clusters information from"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--cells",
        help="Output filename to save cells information",
        type="character", default="cells.tsv"
    )
    parser$add_argument(
        "--umap",
        help="Output filename to save UMAP embeddings",
        type="character", default="umap.tsv"
    )
    parser$add_argument(
        "--clusters",
        help="Output filename to save clusters selected by the --source parameter",
        type="character", default="clusters.tsv"
    )
    args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
    return (args)
}


args <- get_args()

print(paste("Loading Seurat data from", args$rds))
seurat_data <- readRDS(args$rds)

cells_data <- Cells(seurat_data)
print(paste("Loaded", length(cells_data), "cell"))
print(head(cells_data))

umap_data <- Embeddings(seurat_data, reduction="umap")
print(paste("Loaded", length(umap_data), "UMAP embeddings"))
print(head(umap_data))

cluster_data <- seurat_data@meta.data[, args$source]
print(paste("Loaded", length(cluster_data), "clusters"))
print(head(cluster_data))

export_data(cells_data, args$cells, row_names=FALSE, col_names=FALSE)
export_data(umap_data, args$umap, row_names=FALSE, col_names=FALSE)
export_data(cluster_data, args$clusters, row_names=FALSE, col_names=FALSE)