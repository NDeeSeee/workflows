import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("data.table", attach=FALSE)
import("reticulate", attach=FALSE)

export(
    "get_matrix",
    "write_matrix",
    "cb_build",
    "export_cellbrowser"
)


get_matrix <- function(object, slot){
    if (slot == "counts") {
        counts <- SeuratObject::GetAssayData(object=object, slot="counts")
    } else if (slot == "scale.data") {
        counts <- SeuratObject::GetAssayData(object=object, slot="scale.data")
    } else if (slot == "data") {
        counts <- SeuratObject::GetAssayData(object=object)
    } else {
        base::print(base::paste("Slot", slot, "not found. Quit."))
        base::quit(save="no", status=1, runLast=FALSE)
    }
}

write_matrix <- function(inMat, outFname, sliceSize=1000){
    fnames <- c()
    mat <- inMat
    geneCount <- base::nrow(mat)
    startIdx <- 1
    while (startIdx < geneCount){
        endIdx <- min(startIdx + sliceSize - 1, geneCount)
        matSlice <- mat[startIdx:endIdx,]
        denseSlice <- base::as.matrix(matSlice)
        dt <- data.table::data.table(denseSlice)
        dt <- base::cbind(gene=base::rownames(matSlice), dt)
        writeHeader <- FALSE
        if (startIdx == 1) { 
            writeHeader <- TRUE
        }
        sliceFname <- base::paste0("temp", startIdx, ".txt")
        data.table::fwrite(dt, sep="\t", file=sliceFname, quote=FALSE, col.names=writeHeader)
        fnames <- base::append(fnames, sliceFname)
        startIdx <- startIdx + sliceSize
    }
    base::system(base::paste("cat", base::paste(fnames, collapse=" "), "| gzip >", outFname, sep=" "))
    base::unlink(fnames)
}

cb_build <- function(object, slot, rootname, cluster_field, features, meta_fields, meta_fields_names, dot_radius, dot_alpha, hide_labels) {
    if (!base::dir.exists(rootname)) {
        base::dir.create(rootname)
    }
    idents <- SeuratObject::Idents(object)
    meta <- object@meta.data
    cell_order <- base::colnames(object)
    counts <- get_matrix(object, slot)
    genes <- base::rownames(x=object)
    dr <- object@reductions

    gzPath <- base::file.path(rootname, "exprMatrix.tsv.gz")
    if ((((base::ncol(counts)/1000)*(base::nrow(counts)/1000))>2000) && methods::is(counts, "sparseMatrix")){
        write_matrix(counts, gzPath);
    } else {
        mat <- base::as.matrix(counts)
        df <- base::as.data.frame(mat, check.names=FALSE)
        df <- base::data.frame(gene=genes, df, check.names=FALSE)
        z <- base::gzfile(gzPath, "w")
        utils::write.table(x=df, sep="\t", file=z, quote=FALSE, row.names=FALSE)
        base::close(con=z)
    }
    embeddings = names(dr)
    embeddings_conf <- c()
    for (embedding in embeddings) {
        emb <- dr[[embedding]]
        df <- emb@cell.embeddings
        if (base::ncol(df) > 2){
            df <- df[, 1:2]
        }
        base::colnames(df) <- c("x", "y")
        df <- base::data.frame(cellId=base::rownames(df), df, check.names=FALSE)
        utils::write.table(
            df[cell_order,],
            sep="\t",
            file=base::file.path(rootname, base::sprintf("%s.coords.tsv", embedding)),
            quote=FALSE,
            row.names=FALSE
        )
        embeddings_conf <- c(
            embeddings_conf,
            base::sprintf('{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}', embedding)
        )
    }
    df <- base::data.frame(row.names=cell_order, check.names=FALSE)
    enum_fields <- c()
    for (i in 1:length(meta_fields)){
        field <- meta_fields[i]
        name <- meta_fields_names[i]
        if (field %in% base::colnames(meta)){
            df[[name]] <- meta[[field]]
            if (!is.numeric(df[[name]])) {
                enum_fields <- c(enum_fields, name)
            }
        }
    }
    df <- base::data.frame(Cell=base::rownames(df), df, check.names=FALSE)
    utils::write.table(
        base::as.matrix(df[cell_order,]),
        sep="\t",
        file=base::file.path(rootname, "meta.tsv"),
        quote=FALSE,
        row.names=FALSE
    )
    if (!is.null(features)){
        utils::write.table(
            features,
            sep="\t",
            file=base::file.path(rootname, "quickGenes.tsv"),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE
        )
    }
    config <- '
name="cellbrowser"
shortLabel="cellbrowser"
geneLabel="Feature"
exprMatrix="exprMatrix.tsv.gz"
meta="meta.tsv"
radius=%i
alpha=%f
geneIdType="auto"
enumFields=%s
coords=%s'
    enum_string <- base::paste0("[", base::paste(base::paste0('"', enum_fields, '"'), collapse=", "), "]")
    coords_string <- base::paste0("[", base::paste(embeddings_conf, collapse = ",\n"), "]")
    config <- base::sprintf(
        config,
        dot_radius,
        dot_alpha,
        enum_string,
        coords_string
    )
    if (!is.null(features)){
        config <- base::paste(config, 'quickGenesFile="quickGenes.tsv"', sep="\n")
    }
    if (!is.null(cluster_field)){
        config <- base::paste(config, base::paste0('clusterField="', cluster_field,'"'), sep="\n")
        config <- base::paste(config, base::paste0('labelField="', cluster_field,'"'), sep="\n")
    }
    if (!is.null(hide_labels) && hide_labels){
        config <- base::paste(config, 'showLabels=False', sep="\n")
    }
    conf_path = base::file.path(rootname, "cellbrowser.conf")
    base::cat(config, file=conf_path)
    cb_dir <- base::paste(rootname, "html_data", sep="/")
    cb <- reticulate::import(module = "cellbrowser")
    cb$cellbrowser$build(rootname, cb_dir)
}

export_cellbrowser <- function(seurat_data, assay, slot, rootname, dot_radius=3, dot_alpha=0.5, features=NULL, meta_fields=NULL, meta_fields_names=NULL, resolution=NULL, resolution_prefix=NULL){
    base::tryCatch(
        expr = {
            backup_assay <- SeuratObject::DefaultAssay(seurat_data)
            SeuratObject::DefaultAssay(seurat_data) <- assay                 # safety measure
            SeuratObject::Idents(seurat_data) <- "new.ident"                 # safety measure

            if (is.null(meta_fields) || is.null(meta_fields_names)){
                meta_fields <-       c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "S.Score", "G2M.Score",    "Phase", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment",       "nucleosome_signal", "frip", "blacklisted_fraction")
                meta_fields_names <- c("GEX UMIs",   "Genes",        "Mitochondrial %", "Novelty score",            "S score", "G to M score", "Phase", "ATAC UMIs",   "Peaks",         "TSS enrichment score", "Nucleosome signal", "FRiP", "Bl. regions")
            }

            for (i in 1:length(meta_fields)){
                current_field <- meta_fields[i]
                if (current_field %in% base::colnames(seurat_data@meta.data) && length(base::unique(base::as.vector(as.character(seurat_data@meta.data[, current_field])))) > 1){
                    default_cluster_field <- meta_fields_names[i]                 # first field from the meta.data that is not unique for all cells
                    break
                }
            }
            cluster_field <- default_cluster_field                                # default value (in case we have only one identity and no clusters)

            if (length(base::unique(base::as.vector(as.character(seurat_data@meta.data$new.ident)))) > 1){
                cluster_field <- "Identity"
                meta_fields <- base::append(meta_fields, "new.ident", 0)
                meta_fields_names <- base::append(meta_fields_names, "Identity", 0)
            }

            if(seurat_data@meta.data$new.ident != seurat_data@meta.data$condition){
                meta_fields <- base::append(meta_fields, "condition", 1)
                meta_fields_names <- base::append(meta_fields_names, "Condition", 1)
            }

            if(!is.null(resolution) && !is.null(resolution_prefix)){
                meta_fields <- base::append(meta_fields, paste(resolution_prefix, resolution, sep="."))
                meta_fields_names <- base::append(meta_fields_names, paste("Clustering (", resolution, ")", sep=""))
                cluster_field <- paste("Clustering (", resolution, ")", sep="")[1]
            }

            base::tryCatch(
                expr = {
                    custom_fields <- base::grep("^custom_", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
                    custom_fields_names <- base::gsub("custom_", "Custom ", custom_fields)
                    meta_fields <- base::append(meta_fields, custom_fields)
                    meta_fields_names <- base::append(meta_fields_names, custom_fields_names)
                },
                error = function(e){
                    base::print(base::paste("Failed to add custom metadata fields with error -", e))
                }
            )

            base::tryCatch(
                expr = {
                    quartile_fields <- base::grep("^quartile_", base::colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
                    quartile_fields_names <- base::gsub("quartile_", "Quartile ", quartile_fields)
                    meta_fields <- base::append(meta_fields, quartile_fields)
                    meta_fields_names <- base::append(meta_fields_names, quartile_fields_names)
                },
                error = function(e){
                    base::print(base::paste("Failed to add quartile metadata fields with error -", e))
                }
            )

            cb_build(
                seurat_data,
                slot=slot,
                rootname=rootname,
                cluster_field=cluster_field,
                features=features,
                meta_fields=meta_fields,
                meta_fields_names=meta_fields_names,
                dot_radius=dot_radius,
                dot_alpha=dot_alpha,
                hide_labels=base::ifelse(cluster_field==default_cluster_field, TRUE, FALSE)
            )
            base::print(base::paste("Exporting UCSC Cellbrowser data to", rootname, sep=" "))
        },
        error = function(e){
            base::print(base::paste("Failed to export UCSC Cellbrowser data with error -", e))
        },
        finally = {
            SeuratObject::DefaultAssay(seurat_data) <- backup_assay
        }
    )
}