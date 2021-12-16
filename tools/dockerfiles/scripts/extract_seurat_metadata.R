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


export_metadata <- function(data, suffix, args){
    cells_data <- Cells(data)
    print(paste("Saving", length(cells_data), "cell"))
    print(head(cells_data))

    umap_data <- Embeddings(data, reduction="umap")
    print(paste("Saving", length(umap_data), "UMAP embeddings"))
    print(head(umap_data))

    cluster_data <- data@meta.data[, args$source]
    print(paste("Saving", length(cluster_data), "clusters"))
    print(head(cluster_data))

    export_data(
        cells_data,
        location=paste0(args$cells, "_", suffix, ".tsv"),
        row_names=FALSE,
        col_names=FALSE
    )
    export_data(
        umap_data,
        location=paste0(args$umap, "_", suffix, ".tsv"),
        row_names=FALSE,
        col_names=FALSE
    )
    export_data(
        cluster_data,
        location=paste0(args$clusters, "_", suffix, ".tsv"),
        row_names=FALSE,
        col_names=FALSE
    )
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
        "--splitby",
        help=paste(
            "Column name to split data by. Barcode suffixes will be removed.",
            "Output files will have values from the selected column as part of their names.",
            "Default: do not split, export combined file"
        ),
        type="character"
    )
    parser$add_argument(
        "--cells",
        help="Output prefix to save cells information data as TSV file",
        type="character", default="./cells"
    )
    parser$add_argument(
        "--umap",
        help="Output prefix to save UMAP embeddings data as TSV file",
        type="character", default="./umap"
    )
    parser$add_argument(
        "--clusters",
        help=paste(
            "Output prefix to save clusters selected by the --source",
            "parameter data as TSV file"
        ),
        type="character", default="./clusters"
    )
    args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
    return (args)
}


args <- get_args()

print(paste("Loading Seurat data from", args$rds))
seurat_data <- readRDS(args$rds)

if (!is.null(args$splitby)){
    print(paste("Splitting Seurat data by", args$splitby))
    Idents(seurat_data) <- args$splitby
    identities <- unique(as.vector(as.character(Idents(seurat_data))))
    for (i in 1:length(identities)) {
        current_identity <- identities[i]
        print(paste("Subsetting to", current_identity))
        suffix <- gsub("'|\"| |-", "_", current_identity)
        filtered_seurat_data <- subset(seurat_data, idents=current_identity)
        filtered_seurat_data <- RenameCells(
            filtered_seurat_data,
            new.names=sapply( strsplit(Cells(filtered_seurat_data), "-", fixed=TRUE), "[", 1 )
        )
        export_metadata(filtered_seurat_data, suffix, args)
    }
} else {
    print("Exporting no splitted data")
    export_metadata(seurat_data, "combined", args)
}