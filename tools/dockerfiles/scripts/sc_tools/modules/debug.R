import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)

export(
    "print_info"
)


print_info <- function(seurat_data, args) {
    if (!is.null(args$verbose) && args$verbose){
        base::cat("Datasets:\n")
        base::cat(
            base::paste(
                base::unique(base::as.vector(as.character(seurat_data@meta.data$new.ident))),
                collapse="\n"
            )
        )
        base::cat("\n\nConditions:\n")
        base::cat(
            base::paste(
                base::unique(base::as.vector(as.character(seurat_data@meta.data$condition))),
                collapse="\n"
            )
        )
        base::cat("\n\nAssays:\n")
        base::print(seurat_data@assays)
        base::cat("Default assay:\n")
        base::cat(SeuratObject::DefaultAssay(seurat_data))
        base::cat("\n\nReductions:\n")
        base::print(seurat_data@reductions)
        base::cat("Graphs:\n")
        base::cat(
            base::paste(
                base::unique(base::as.vector(as.character(SeuratObject::Graphs(seurat_data)))),
                collapse="\n"
            )
        )
        base::cat("\n\nNeighbors:\n")
        base::cat(
            base::paste(
                base::unique(base::as.vector(as.character(SeuratObject::Neighbors(seurat_data)))),
                collapse="\n"
            )
        )
        base::cat("\n\nMiscellaneous:\n")
        base::print(SeuratObject::Misc(seurat_data))
        if ("ATAC" %in% names(seurat_data@assays)) {
            base::cat("Fragments from ATAC assay:\n")
            fragments <- Signac::Fragments(seurat_data[["ATAC"]])
            if (length(fragments) > 0){
                for (i in 1:length(fragments)){
                    fragment <- fragments[[i]]
                    base::cat(
                        base::paste(
                            "Fragment", i, "includes",
                            length(methods::slot(fragment, "cells")), "cells",
                            "loaded from", methods::slot(fragment, "path"),
                            "with hash", methods::slot(fragment, "hash")[1],       # hash 1 is for Fragment file, hash 2 is for index
                            "\n"
                        )
                    )
                }
            }
        }
        base::cat("\nMetadata:\n")
        base::print(utils::head(seurat_data@meta.data))
        base::cat("\nMemory usage:\n")
        base::print(base::gc())
    } else {
        base::gc(verbose=FALSE)                                                    # silently cleans gc
    }
}