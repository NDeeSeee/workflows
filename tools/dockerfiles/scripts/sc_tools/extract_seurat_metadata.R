#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(argparse))
suppressMessages(library(data.table))


setup_parallelization <- function (args) {
    invisible(capture.output(plan("multiprocess", workers=args$cpus)))
    invisible(capture.output(plan()))
    invisible(capture.output(setDTthreads(args$cpus)))
    options(future.globals.maxSize = args$memory * 1024^3)               # convert to bytes
}


debug_seurat_data <- function(seurat_data, args) {
    if (args$verbose){
        print("Assays")
        print(seurat_data@assays)
        print("Reductions")
        print(seurat_data@reductions)
        print("Active assay")
        print(DefaultAssay(seurat_data))
        print("Metadata")
        print(head(seurat_data@meta.data))
    }
}


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


export_metadata <- function(data, suffix, args){
    cluster_data <- cbind(Cells(data), as.vector(data@meta.data[, args$source]))  # need to make sure we don't take levels instead of the values. Otherwise the order is not correct
    export_data(
        cluster_data,
        location=paste(args$output, suffix, "clusters.tsv", sep="_"),
        row_names=FALSE,
        col_names=FALSE
    )
    for (i in 1:length(args$reductions)) {
        current_reduction <- args$reductions[i]
        export_data(
            Embeddings(data, reduction=current_reduction),                                      # if raises exception, it will be caught within export_data
            location=paste(args$output, suffix, paste0(current_reduction, ".tsv"), sep="_"),
            row_names=FALSE,
            col_names=FALSE
        )
    }
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
            print(paste("Failed to export data to", location, "due to", e, sep=" "))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description="Extracts cells, UMAP, and clusters from the Seurat RDS file")
    parser$add_argument(
        "--rds",
        help=paste(
            "Path to the RDS file to load Seurat object from. RDS file can be produced by run_seurat.R,",
            "assign_cell_types.R, or run_seurat_wnn.R scripts."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--source",
        help=paste(
            "Column from metadata to extract clusters information from."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from metadata to split output files by. Barcode suffixes that define",
            "the origin of the cells will be removed. Values from this column will be",
            "appended to the names of output files.",
            "Default: do not split, export combined file"
        ),
        type="character"
    )
    parser$add_argument(
        "--reductions",
        help=paste(
            "List of reductions to extract metadata from. Each reduction",
            "name will be appended to the output files names.",
            "Default: umap, pca, atac_lsi, rnaumap, atacumap, wnnumap"
        ),
        type="character",
        default=c("umap", "pca", "atac_lsi", "rnaumap", "atacumap", "wnnumap"),
        nargs="*"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./metadata",
        type="character", default="./metadata"
    )
    parser$add_argument(
        "--verbose",
        help="Print debug information. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--cpus",
        help="Number of cores/cpus to use. Default: 1",
        type="integer", default=1
    )
    parser$add_argument(
        "--memory",
        help=paste(
            "Maximum memory in GB allowed to be shared between the workers",
            "when using multiple --cpus.",
            "Default: 32"
        ),
        type="integer", default=32
    )
    args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
    return (args)
}


args <- get_args()
print(args)

print(paste("Setting parallelization parameters to", args$cpus, "cores, and", args$memory, "GB of memory"))
setup_parallelization(args)

print(paste("Loading Seurat data from", args$rds))
seurat_data <- readRDS(args$rds)
debug_seurat_data(seurat_data, args)

if (!is.null(args$splitby)){
    print(paste("Exporting results splitted by", args$splitby))
    Idents(seurat_data) <- args$splitby
    identities <- unique(as.vector(as.character(Idents(seurat_data))))
    for (i in 1:length(identities)) {
        current_identity <- identities[i]
        print(paste("Processing", current_identity))
        suffix <- gsub("'|\"| |-", "_", current_identity)                       # safety measure
        filtered_seurat_data <- subset(seurat_data, idents=current_identity)
        filtered_seurat_data <- RenameCells(
            filtered_seurat_data,
            new.names=sapply( strsplit(Cells(filtered_seurat_data), "-", fixed=TRUE), "[", 1 )
        )
        export_metadata(filtered_seurat_data, suffix, args)
    }
} else {
    print("Exporting combined results")
    export_metadata(seurat_data, "combined", args)
}