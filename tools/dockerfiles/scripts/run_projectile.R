#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB should be good by default

suppressMessages(library(future))
suppressMessages(library(argparse))
suppressMessages(library(ProjecTILs))


set_threads <- function (threads) {
    invisible(capture.output(plan("multiprocess", workers=threads)))
    invisible(capture.output(plan()))
    invisible(capture.output(setDTthreads(threads)))
}


export_projection_plot <- function(reference_data, projected_data, rootname, plot_title, labels_col="functional.cluster", linesize=1, pointsize=1, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {

            plot <- plot.projection(
                        ref=reference_data,
                        query=projected_data,
                        labels.col=labels_col,
                        linesize=linesize,
                        pointsize=pointsize
                    ) +
                    theme_gray() +
                    ggtitle(plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export projection plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export projection plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_composition_plot <- function(reference_data, projected_data, rootname, plot_title, labels_col="functional.cluster", metric="Percent", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {

            plot <- plot.statepred.composition(
                        ref=reference_data,
                        query=projected_data,
                        labels.col=labels_col,
                        metric=metric
                    ) +
                    theme_gray() +
                    ggtitle(plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export composition plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export composition plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_radar_plot <- function(reference_data, projected_data, rootname, labels_col="functional.cluster", min_cells=30, genes4radar=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {

            plot <- plot.states.radar(
                        ref=reference_data,
                        query=projected_data,
                        labels.col=labels_col,
                        genes4radar=genes4radar,
                        min.cells=min_cells,
                        return=TRUE                    # need it to return as plot
                    )

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(plot(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(plot(plot))
                dev.off()
            }

            print(paste("Export radar plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to radar composition plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


run_projection <- function(reference_data, query_data, suffix, args){
    projected_data <- make.projection(
        query_data,
        ref=reference_data,
        filter.cells=!args$nofilter,
        query.assay="RNA",
        skip.normalize=TRUE       # our Seurat data is already normalized
    )
    export_projection_plot(
        reference_data=reference_data,
        projected_data=projected_data,
        rootname=paste(args$output, "projection", suffix, sep="_"),
        plot_title=paste("Projection of group", suffix, "into reference UMAP")
    )
    projected_data <- cellstate.predict(ref=reference_data, query=projected_data)
    print(table(projected_data$functional.cluster))
    export_composition_plot(
        reference_data=reference_data,
        projected_data=projected_data,
        rootname=paste(args$output, "composition", suffix, sep="_"),
        plot_title=paste("Cell state composition of group", suffix)
    )
    export_radar_plot(
        reference_data=reference_data,
        projected_data=projected_data,
        genes4radar=args$features,
        rootname=paste(args$output, "radar", suffix, sep="_")
    )
}


get_args <- function(){
    parser <- ArgumentParser(description="Projects query Seurat data into reference single-cell atlas")
    parser$add_argument(
        "--query",
        help=paste(
            "Path to the RDS file to load query Seurat object from",
            "(RNA assay should be already normalized). RDS file can",
            "be produced by run_seurat.R or assign_cell_types.R scripts."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--reference",
        help=paste(
            "Path to the RDS file to load reference Seurat object from.",
            "Default: referenceTIL downloaded on the fly"
        ),
        type="character", default="referenceTIL"
    )
    parser$add_argument(
        "--splitby",
        help=paste(
            "Column from query Seurat object metadata to split projection by.",
            "Default: do not split, process entire query dataset at once"
        ),
        type="character"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./projectils",
        type="character", default="./projectils"
    )
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--features",
        help="Features of interest to highlight. Default: Foxp3, Cd4, Cd8a, Tcf7, Ccr7, Gzmb, Gzmk, Pdcd1, Havcr2, Tox, Mki67, Il2ra",
        type="character", nargs="*", default=c("Foxp3", "Cd4", "Cd8a", "Tcf7", "Ccr7", "Gzmb", "Gzmk", "Pdcd1", "Havcr2", "Tox", "Mki67", "Il2ra")
    )
    parser$add_argument(
        "--nofilter",
        help="Skip pre-filtering T cells using scGate. Default: pre-filter T-cells",
        action="store_true"
    )
    parser$add_argument(
        "--threads",
        help="Threads. Default: 1",
        type="integer", default=1
    )
    args <- parser$parse_args(commandArgs(trailingOnly=TRUE))
    return (args)
}


args <- get_args()
print("Used parameters")
print(args)
print(paste("Setting parallelizations threads to", args$threads))
set_threads(args$threads)                                                  # not sure if it actually speeds up make.projection

print(paste("Loading reference data from", args$ref))
reference_data <- load.reference.map(ref=args$ref)

print(paste("Loading query data from", args$query))
query_data <- readRDS(args$query)
DefaultAssay(query_data) <- "RNA"                                          # safety measure as we need to work with RNA assay

if (!is.null(args$splitby)){
    print(paste("Exporting results splitted by", args$splitby))
    Idents(query_data) <- args$splitby
    groups <- unique(as.vector(as.character(Idents(query_data))))
    for (i in 1:length(groups)) {
        current_group <- groups[i]
        print(paste("Processing", current_group))
        suffix <- gsub("'|\"| |-", "_", current_group)                     # safety measure
        filtered_query_data <- subset(query_data, idents=current_group)
        tryCatch(
            expr = {
                run_projection(reference_data, filtered_query_data, suffix, args)
            },
            error = function(e){
                print(paste("Failed to run projection for", suffix))
            }
        )
    }
} else {
    print("Exporting combined results")
    run_projection(reference_data, filtered_query_data, "combined", args)
}