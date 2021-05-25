#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)

suppressMessages(library(tidyverse))
suppressMessages(library(argparse))


export_data <- function(data, location, row_names=FALSE, col_names=TRUE, quote=FALSE){
    write.table(
        data,
        file=location,
        sep="\t",
        row.names=row_names,
        col.names=col_names,
        quote=quote
    )
}


export_garnett_formatted_data <- function(data, location){
    for(i in 1:nrow(data)) {
        row <- data[i,]
        cat(paste(">", row[1], "\n"), file=location, append=TRUE)
        cat(paste("expressed:", row[2], "\n"), file=location, append=TRUE)
    }
}


load_cell_markers_data <- function (location, species, gene_sep=",") {
    raw_cell_markers_data <- read.table(
        location,
        sep="\t",
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    exp_columns <- paste0(
        "exp_",
        seq_len(
            max(length(strsplit(raw_cell_markers_data$geneSymbol, gene_sep)))
        )
    )
    formatted_cell_markers_data <- raw_cell_markers_data %>%
                                   filter(.$speciesType==species) %>%
                                   separate("geneSymbol", all_of(exp_columns), gene_sep) %>%
                                   gather("not_used", "Marker", all_of(exp_columns), na.rm=TRUE) %>%
                                   select("cellName", "Marker") %>%
                                   rename("Type"="cellName") %>%
                                   mutate("Marker"=str_trim(.$Marker)) %>%
                                   filter(.$Marker!="NA") %>%
                                   distinct() %>%
                                   group_by(Type) %>%
                                   mutate(Expressed=paste0(Marker, collapse=", ")) %>%
                                   select("Type", "Expressed") %>%
                                   distinct()
    return (formatted_cell_markers_data)
}


get_args <- function(){
    parser <- ArgumentParser(description='Converts single cell markers from CellMarker to Garnett format')
    parser$add_argument("--raw",     help="Path to the raw file from http://biocc.hrbmu.edu.cn/CellMarker/download.jsp", type="character", required="True")
    parser$add_argument("--species", help="Select species to save", type="character", choices=c("Human", "Mouse"), required="True")
    parser$add_argument("--output",  help="Output prefix. Default: ./garnett", type="character", default="./garnett")
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}


args <- get_args()

print("Loading cell markers data")
cell_markers_data <- load_cell_markers_data(args$raw, args$species)

print("Exporting cell markers data")
export_data(cell_markers_data, paste(args$output, tolower(args$species), "markers.tsv", sep="_"))
export_garnett_formatted_data(cell_markers_data, paste(args$output, tolower(args$species), "formatted_markers.tsv", sep="_"))




