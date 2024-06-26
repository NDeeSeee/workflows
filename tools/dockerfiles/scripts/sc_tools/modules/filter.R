import("Seurat")
import("Signac")

export(
    "remove_doublets",
    "apply_rna_qc_filters",
    "apply_atac_qc_filters",
    "apply_peak_qc_filters",
    "collapse_fragments_list"
)


collapse_fragments_list <- function(seurat_data){
    if ("ATAC" %in% names(seurat_data@assays)) {
        fragments_list <- Signac::Fragments(seurat_data[["ATAC"]])
        if (length(fragments_list) > 1){
            base::print("Collapsing repetitive fragment objects")
            collapsed_fragments <- fragments_list[[1]]                                       # will use the 1 Fragment assuming that all of them reference the same file
            collapsed_hash <- methods::slot(collapsed_fragments, "hash")[1]                  # if hash is not the same for all Fragment objects, quit
            collected_cells <- c()
            for (i in 1:length(fragments_list)){
                current_fragments <- fragments_list[[i]]
                current_hash <- methods::slot(current_fragments, "hash")[1]
                current_path <- methods::slot(current_fragments, "path")
                current_cells <- Signac::Cells(current_fragments)
                base::print(
                    base::paste(
                        "Fragments object", i, "includes",
                        length(current_cells), "cells loaded from",
                        current_path, "with hash", current_hash
                    )
                )
                collected_cells <- base::unique(c(collected_cells, current_cells))           # use unique as a safety measure
                if (collapsed_hash != current_hash){
                    base::print("Attempting to collapse different Fragment objects. Quit.")
                    base::quit(save="no", status=1, runLast=FALSE)
                }
                base::rm(current_fragments, current_hash, current_cells)
            }
            names(collected_cells) <- collected_cells                                        # need to be named vector
            Signac::Cells(collapsed_fragments) <- collected_cells
            seurat_data[["ATAC"]] <- SeuratObject::SetAssayData(
                seurat_data[["ATAC"]],
                slot="fragments",
                new.data=collapsed_fragments
            )
            base::rm(collapsed_fragments, collapsed_hash, collected_cells)
        }
        base::rm(fragments_list)
    }
    return (seurat_data)
}

# apply_cell_filters <- function(seurat_data, barcodes_data) {
#     filtered_seurat_data <- base::subset(seurat_data, cells=barcodes_data)
#     return (filtered_seurat_data)
# }

remove_doublets <- function(seurat_data, what_to_remove) {
    base::print(base::paste("Cells before filtering:", length(SeuratObject::Cells(seurat_data))))
    if (what_to_remove == "union"){
        base::print("Removing doublets identified in either of RNA or ATAC assays")
        seurat_data <- base::subset(
            seurat_data,
            subset=(rna_doublets == "singlet") & (atac_doublets == "singlet")
        )
    } else if (what_to_remove == "onlyrna"){
        base::print("Removing only RNA doublets")
        seurat_data <- base::subset(
            seurat_data,
            subset=(rna_doublets == "singlet")
        )
    } else if (what_to_remove == "onlyatac"){
        base::print("Removing only ATAC doublets")
        seurat_data <- base::subset(
            seurat_data,
            subset=(atac_doublets == "singlet")
        )
    } else if (what_to_remove == "intersect"){
        base::print("Removing doublets identified as common for RNA and ATAC assays")
        seurat_data <- base::subset(
            seurat_data,
            subset=(rna_doublets == "singlet") | (atac_doublets == "singlet")
        )
    }
    base::print(base::paste("Cells after filtering:", length(SeuratObject::Cells(seurat_data))))
    return (seurat_data)
}

apply_rna_qc_filters <- function(seurat_data, args) {
    base::print(base::paste("Cells before filtering:", length(SeuratObject::Cells(seurat_data))))
    merged_seurat_data <- NULL
    SeuratObject::Idents(seurat_data) <- "new.ident"                                                 # safety measure
    sorted_identities <- base::sort(                                                                 # alphabetically sorted identities A -> Z
        base::unique(
            base::as.vector(as.character(SeuratObject::Idents(seurat_data)))
        )
    )
    for (i in 1:length(sorted_identities)){
        identity <- sorted_identities[i]

        mingenes <- args$mingenes[i]
        maxgenes <- args$maxgenes[i]
        minumis <- args$minumis[i]
        minnovelty <- args$minnovelty[i]

        base::print(base::paste("Filtering", identity))
        base::print(base::paste(" ", mingenes, "<= Genes per cell <=", maxgenes))
        base::print(base::paste(" ", "RNA reads per cell >=", minumis))
        base::print(base::paste(" ", "Novelty score >=", minnovelty))
        base::print(base::paste(" ", "Percentage of RNA reads mapped to mitochondrial genes <=", args$maxmt))

        filtered_seurat_data <- base::subset(
            seurat_data,
            idents=identity,
            subset=(nFeature_RNA >= mingenes) &
                   (nFeature_RNA <= maxgenes) &
                   (nCount_RNA >= minumis) &
                   (log10_gene_per_log10_umi >= minnovelty) &
                   (mito_percentage <= args$maxmt)
        )

        if (is.null(merged_seurat_data)){
            merged_seurat_data <- filtered_seurat_data
        } else {
            merged_seurat_data <- base::merge(merged_seurat_data, y=filtered_seurat_data)
        }
        base::rm(filtered_seurat_data)                                                          # remove unused data
    }
    base::print(base::paste("Cells after filtering:", length(SeuratObject::Cells(merged_seurat_data))))
    merged_seurat_data <- collapse_fragments_list(merged_seurat_data)
    base::gc(verbose=FALSE)
    return (merged_seurat_data)
}

apply_atac_qc_filters <- function(seurat_data, args) {
    base::print(base::paste("Cells before filtering:", length(SeuratObject::Cells(seurat_data))))
    merged_seurat_data <- NULL
    SeuratObject::Idents(seurat_data) <- "new.ident"                                                 # safety measure
    sorted_identities <- sort(unique(as.vector(as.character(Idents(seurat_data)))))                  # alphabetically sorted identities A -> Z
    for (i in 1:length(sorted_identities)){
        identity <- sorted_identities[i]

        minfragments <- args$minfragments[i]
        maxnuclsignal <- args$maxnuclsignal[i]
        mintssenrich <- args$mintssenrich[i]

        base::print(base::paste("Filtering", identity))
        base::print(base::paste(" ", "ATAC fragments in peaks per cell >=", minfragments))
        base::print(base::paste(" ", "Nucleosome signal <=", maxnuclsignal))
        base::print(base::paste(" ", "TSS enrichment score >=", mintssenrich))

        filtered_seurat_data <- base::subset(
            seurat_data,
            idents=identity,
            subset=(nCount_ATAC >= minfragments) &
                   (nucleosome_signal <= maxnuclsignal) &
                   (TSS.enrichment >= mintssenrich)
        )

        if (is.null(merged_seurat_data)){
            merged_seurat_data <- filtered_seurat_data
        } else {
            merged_seurat_data <- base::merge(merged_seurat_data, y=filtered_seurat_data)
        }
        base::rm(filtered_seurat_data)                                                          # remove unused data
    }
    base::print(base::paste("Cells after filtering:", length(SeuratObject::Cells(merged_seurat_data))))
    merged_seurat_data <- collapse_fragments_list(merged_seurat_data)
    base::gc(verbose=FALSE)
    return (merged_seurat_data)
}

apply_peak_qc_filters <- function(seurat_data, args) {
    base::print(base::paste("Cells before filtering:", length(SeuratObject::Cells(seurat_data))))
    merged_seurat_data <- NULL
    SeuratObject::Idents(seurat_data) <- "new.ident"                                                 # safety measure
    sorted_identities <- sort(unique(as.vector(as.character(Idents(seurat_data)))))                  # alphabetically sorted identities A -> Z
    for (i in 1:length(sorted_identities)){
        identity <- sorted_identities[i]

        maxblacklist <- args$maxblacklist[i]

        base::print(base::paste("Filtering", identity))
        base::print(base::paste(" ", "FRiP >=", args$minfrip))
        base::print(base::paste(" ", "Ratio of reads in genomic blacklist regions <=", maxblacklist))

        filtered_seurat_data <- base::subset(
            seurat_data,
            idents=identity,
            subset=(frip >= args$minfrip) & (blacklist_fraction <= maxblacklist)
        )
        if (is.null(merged_seurat_data)){
            merged_seurat_data <- filtered_seurat_data
        } else {
            merged_seurat_data <- base::merge(merged_seurat_data, y=filtered_seurat_data)
        }
        base::rm(filtered_seurat_data)                                                       # remove unused data
    }
    base::print(base::paste("Cells after filtering:", length(SeuratObject::Cells(merged_seurat_data))))
    merged_seurat_data <- collapse_fragments_list(merged_seurat_data)
    base::gc(verbose=FALSE)
    return (merged_seurat_data)
}