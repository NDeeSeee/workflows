#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(argparse))
suppressMessages(library(tidyverse))
suppressMessages(library(UpSetR))

get_args <- function(){
    parser <- ArgumentParser(description="Souporcell Cluster Overlap")
    parser$add_argument(
        "--clusters",
        help="Path to the Souporcell generated cluster files",
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--aliases",
        help=paste(
            "Aliases for cluster files. Order and number of the",
            "provided values should correspond to --clusters."
        ),
        type="character", required="True", nargs="+"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./overlap",
        type="character", default="./overlap"
    )
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}

args <- get_args()

print("Input parameters")
print(args)

collected_list <- list()

for (i in 1:length(args$clusters)){
    current_location <- args$clusters[i]
    current_alias <- args$aliases[i]
    print(paste("Loading from", current_location, "as", current_alias))
    cluster_data <- utils::read.table(
                        current_location,
                        sep="\t",
                        header=TRUE,
                        check.names=FALSE,
                        stringsAsFactors=FALSE
                    ) %>%
                    dplyr::filter(.$status == "singlet") %>%
                    dplyr::select(c("barcode", "assignment")) %>%
                    dplyr::group_by(assignment)
    split_data <- dplyr::group_split(cluster_data)
    group_names <- dplyr::group_keys(cluster_data) %>% pull("assignment")
    for (j in 1:length(group_names)){
        current_group_name <- group_names[j]
        collected_list[[paste(current_alias, as.character(current_group_name), sep="_")]] <- split_data[j][[1]] %>% pull("barcode")
    }
}

grDevices::png(
    file=paste0(args$output, ".png"),
    width=1200,
    height=400,
    res=100
)
print(upset(fromList(collected_list), order.by="freq"))
grDevices::dev.off()