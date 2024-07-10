#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(modules))
suppressMessages(library(argparse))

HERE <- (function() {return (dirname(sub("--file=", "", commandArgs(trailingOnly=FALSE)[grep("--file=", commandArgs(trailingOnly=FALSE))])))})()
suppressMessages(io <- modules::use(file.path(HERE, "modules/io.R")))
suppressMessages(prod <- modules::use(file.path(HERE, "modules/prod.R")))


get_args <- function(){
    parser <- ArgumentParser(description="Creates ATAC Search file required by UCSC Cellbrowser")
    parser$add_argument(
        "--annotations",
        help="Path to the genome annotation file in GTF format",
        type="character", required="True"
    )
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    print(args)
    return (args)
}


args <- get_args()

print(paste("Loading genome annotation data from", args$annotations))
annotation_data <- base::as.data.frame(io$load_annotation_data(args$annotations)) %>%
                   dplyr::filter(.$type=="transcript") %>%
                   dplyr::group_by(seqnames, strand, gene_id) %>%                          # to account for a replicates in the gene names
                   dplyr::arrange(desc(width), .by_group=TRUE) %>%
                   dplyr::slice_head(n=1) %>%                                              # takes the coordinates of the longest transcript within each group
                   dplyr::ungroup() %>%
                   dplyr::mutate("score"=0) %>%                                            # add empty score column to corresond to BED format
                   dplyr::select(
                       c("seqnames", "start", "end", "gene_id", "score", "strand")
                   )
gene_symbols <- annotation_data %>%
                dplyr::select(c("gene_id")) %>%
                dplyr::mutate("gene_id_copy"=gene_id)

cb_dir <- base::path.expand("~/cellbrowserData/genes")                                     # should be saved in the user's home directory
if (!base::dir.exists(cb_dir)) {
    base::dir.create(cb_dir, recursive=TRUE)
}
utils::write.table(
    annotation_data,
    base::gzfile(base::file.path(cb_dir, "genome.current.bed.gz")),
    sep="\t",
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE
)
utils::write.table(
    gene_symbols,
    base::gzfile(base::file.path(cb_dir, "current.symbols.tsv.gz")),
    sep="\t",
    row.names=FALSE,
    col.names=FALSE,
    quote=FALSE
)
exit_code <- sys::exec_wait(
    cmd="cbGenes",
    args=c(
        "json",
        "genome",
        "current",
        base::path.expand(base::file.path(cb_dir, "genome.current.json"))
    )
)
if (exit_code != 0){
    base::print(
        base::paste0(
            "Failed to create genome.current.json ",
            "with exit code ", exit_code, ". Exiting."
        )
    )
    base::quit(save="no", status=1, runLast=FALSE)
}