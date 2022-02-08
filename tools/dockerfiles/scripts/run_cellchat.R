#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB should be good by default
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(argparse))
suppressMessages(library(CellChat))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))


set_threads <- function (threads) {
    invisible(capture.output(plan("multiprocess", workers=threads)))
    invisible(capture.output(plan()))
    invisible(capture.output(setDTthreads(threads)))
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


run_cellchat <- function (data, args, assay="RNA", ctypes=NULL) {
    print("Running CellChat for selected cell types")
    backup_assay <- DefaultAssay(data)
    DefaultAssay(data) <- assay

    Idents(data) <-  args$source
    default_ctypes <- unique(as.vector(as.character(Idents(data))))
    if (is.null(ctypes)){
        ctypes <- default_ctypes
    }
    data_selected <- subset(data, idents=ctypes)
    Idents(data) <- "new.ident"

    print("Creating cellchat object")
    cellchat_data <- createCellChat(data_selected, group.by=args$source)
    print("Setting the ligand-receptor interaction database")
    CellChatDB <- CellChatDB.mouse
    print("Subsetting to only Secreted Signaling")
    CellChatDB.use <- subsetDB(CellChatDB, search="Secreted Signaling")
    cellchat_data@DB <- CellChatDB.use
    cellchat_data <- subsetData(cellchat_data)
    cellchat_data <- identifyOverExpressedGenes(cellchat_data)
    cellchat_data <- identifyOverExpressedInteractions(cellchat_data)
    cellchat_data <- projectData(cellchat_data, PPI.mouse)
    cellchat_data <- computeCommunProb(cellchat_data)
    cellchat_data <- filterCommunication(cellchat_data, min.cells=10)

    cellchat_data <- computeCommunProbPathway(cellchat_data)
    cellchat_data <- aggregateNet(cellchat_data)
    groupSize <- as.numeric(table(cellchat_data@idents))

    par(mfrow = c(1,2), xpd=TRUE)
    pdf(file=paste(paste(args$output, "all_interactions", sep="_"), ".pdf", sep=""))
    netVisual_circle(cellchat_data@net$count, vertex.weight=groupSize, weight.scale=T, label.edge=F, title.name="Number of interactions")
    netVisual_circle(cellchat_data@net$weight, vertex.weight=groupSize, weight.scale=T, label.edge=F, title.name="Interaction weights")
    mat <- cellchat_data@net$weight
    par(mfrow = c(3,4), xpd=TRUE)
    for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight=groupSize, weight.scale=T, edge.weight.max=max(mat), title.name=rownames(mat)[i])
    }
    dev.off()

    # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    # Compute the network centrality scores
    pdf(file=paste(paste(args$output, "network", sep="_"), ".pdf", sep=""))
    cellchat_data <- netAnalysis_computeCentrality(cellchat_data, slot.name="netP")
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat_data, pattern="outgoing")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat_data, pattern="incoming")
    print(ht1)
    print(ht2)
    dev.off()

    df_net_pathways <- subsetCommunication(cellchat_data, slot.name="netP")
    pathways_to_show <- unique(df_net_pathways$pathway_name)
    print("Pathways to show")
    print(pathways_to_show)

    for (i in 1:length(pathways_to_show)) {
        pathway <- pathways_to_show[i]
        pdf(file=paste(paste(args$output, "pathway", pathway, sep="_"), ".pdf", sep=""))
        # Chord diagram
        par(mfrow=c(1,1))
        netVisual_aggregate(cellchat_data, signaling=pathway, layout="chord")
        # Circle plot
        par(mfrow=c(1,1))
        netVisual_aggregate(cellchat_data, signaling=pathway, layout="circle")
        # Computed centrality scores using heatmap
        print(netAnalysis_signalingRole_network(cellchat_data, signaling=pathway, width=8, height=2.5, font.size=10))
        # Contribution of each ligand-receptor pair
        plot <- netAnalysis_contribution(cellchat_data, signaling=pathway)
        print(plot)
        plot <- plotGeneExpression(cellchat_data, signaling=pathway)
        print(plot)
        pairs_lr <- extractEnrichedLR(cellchat_data, signaling=pathway, geneLR.return=FALSE)
        print(pairs_lr)
        print(length(pairs_lr))
        print(rownames(pairs_lr))
        for (i in 1:length(rownames(pairs_lr))) {
            # Circle plot for each ligand-receptor pair
            lr_to_show <- pairs_lr[i,]
            print(lr_to_show)
            print(netVisual_individual(cellchat_data, signaling=pathway, pairLR.use=lr_to_show, layout="circle"))
        }
        dev.off()
    }
    DefaultAssay(data) <- backup_assay
}


get_filtered_data <- function(seurat_data, args){
    if (!is.null(args$groupby) && !is.null(args$select)){
        print(
            paste0(
                "Include only '", paste(args$select, collapse=" "),
                "' values from the '", args$groupby, "' metadata column"
            )
        )
        print(paste("Cells before filtering", nrow(seurat_data@meta.data)))
        Idents(seurat_data) <- args$groupby
        seurat_data <- subset(seurat_data, idents=args$select)
        Idents(seurat_data) <- "new.ident"
        print(paste("Cells after filtering", nrow(seurat_data@meta.data)))
    }
    return (seurat_data)
}


get_args <- function(){
    parser <- ArgumentParser(description="Inference and analysis of cell-cell communication")
    parser$add_argument(
        "--rds",
        help=paste(
            "Path to the RDS file to load Seurat object from. RDS file produced by run_seurat.R script"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--source",
        help=paste(
            "Column name to select clusters for running CellChat"
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--groupby",
        help=paste(
            "Column from the Seurat object metadata to group cells for optional",
            "subsetting with --select parameter (for example, subset to the",
            "specific dataset or condition).",
            "Default: do not subset"
        ),
        type="character"
    )
    parser$add_argument(
        "--select",
        help=paste(
            "Value(s) from the column set with --groupby to optionally subset cells",
            "before running CellChat analysis.",
            "Default: do not subset, use all cells."
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--verbose",
        help="Print debug information. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./cellchat",
        type="character", default="./cellchat"
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

print("Used parameters")
print(args)
print(paste("Setting parallelizations threads to", args$threads))
set_threads(args$threads)

print("Loading data")
seurat_data <- readRDS(args$rds)
debug_seurat_data(seurat_data, args)
seurat_data <- get_filtered_data(seurat_data, args)
debug_seurat_data(seurat_data, args)

run_cellchat(
    seurat_data,
    args
)