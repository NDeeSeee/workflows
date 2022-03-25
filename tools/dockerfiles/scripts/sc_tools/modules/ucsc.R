import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("data.table", attach=FALSE)
import("reticulate", attach=FALSE)

export(
    "cellbrowser_build"
)


find_matrix <- function(object, matrix.slot){
    if (matrix.slot == "counts") {
        counts <- SeuratObject::GetAssayData(object=object, slot="counts")
    } else if (matrix.slot == "scale.data") {
        counts <- SeuratObject::GetAssayData(object=object, slot="scale.data")
    } else if (matrix.slot == "data") {
        counts <- SeuratObject::GetAssayData(object=object)
    } else {
        error("matrix.slot can only be one of: counts, scale.data, data")
    }
}


write_sparse_matrix <- function(inMat, outFname, sliceSize=1000){
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


cellbrowser_build <- function(
    object,
    matrix.slot,
    dir,
    cluster.field,
    features,
    meta.fields,
    meta.fields.names
) {
    if (!base::dir.exists(dir)) {
        base::dir.create(dir)
    }
    idents <- SeuratObject::Idents(object)
    meta <- object@meta.data
    cellOrder <- base::colnames(object)
    counts <- find_matrix(object, matrix.slot)
    genes <- base::rownames(x=object)
    dr <- object@reductions

    gzPath <- base::file.path(dir, "exprMatrix.tsv.gz")
    if ((((base::ncol(counts)/1000)*(base::nrow(counts)/1000))>2000) && methods::is(counts, 'sparseMatrix')){
        write_sparse_matrix(counts, gzPath);
    } else {
        mat <- base::as.matrix(counts)
        df <- base::as.data.frame(mat, check.names=FALSE)
        df <- base::data.frame(gene=genes, df, check.names=FALSE)
        z <- base::gzfile(gzPath, "w")
        utils::write.table(x=df, sep="\t", file=z, quote=FALSE, row.names=FALSE)
        base::close(con=z)
    }
    embeddings = names(dr)
    embeddings.conf <- c()
    for (embedding in embeddings) {
        emb <- dr[[embedding]]
        df <- emb@cell.embeddings
        if (base::ncol(df) > 2){
            df <- df[, 1:2]
        }
        base::colnames(df) <- c("x", "y")
        df <- base::data.frame(cellId=base::rownames(df), df, check.names=FALSE)
        utils::write.table(
            df[cellOrder,],
            sep="\t",
            file=base::file.path(dir, base::sprintf("%s.coords.tsv", embedding)),
            quote=FALSE,
            row.names=FALSE
        )
        embeddings.conf <- c(
            embeddings.conf,
            base::sprintf('{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}', embedding)
        )
    }
    df <- base::data.frame(row.names=cellOrder, check.names=FALSE)
    enum.fields <- c()
    for (i in 1:length(meta.fields)){
        field <- meta.fields[i]
        name <- meta.fields.names[i]
        if (field %in% base::colnames(meta)){
            df[[name]] <- meta[[field]]
            if (!is.numeric(df[[name]])) {
                enum.fields <- c(enum.fields, name)
            }
        }
    }
    df <- base::data.frame(Cell=base::rownames(df), df, check.names=FALSE)
    utils::write.table(
        base::as.matrix(df[cellOrder,]),
        sep="\t",
        file=base::file.path(dir, "meta.tsv"),
        quote=FALSE,
        row.names=FALSE
    )
    if (!is.null(features)){
        utils::write.table(
            features,
            sep="\t",
            file=base::file.path(dir, "quickGenes.tsv"),
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE
        )
    }
    config <- '
name="cellbrowser"
shortLabel="cellbrowser"
geneLabel="Gene"
exprMatrix="exprMatrix.tsv.gz"
meta="meta.tsv"
radius=3
alpha=0.5
geneIdType="auto"
clusterField="%s"
labelField="%s"
enumFields=%s
coords=%s'
    enum.string <- base::paste0("[", base::paste(base::paste0('"', enum.fields, '"'), collapse=", "), "]")
    coords.string <- base::paste0("[", base::paste(embeddings.conf, collapse = ",\n"), "]")
    config <- base::sprintf(
        config,
        cluster.field,
        cluster.field,
        enum.string,
        coords.string
    )
    if (!is.null(features)){
        config <- base::paste(config, 'quickGenesFile="quickGenes.tsv"', sep="\n")
    }
    confPath = base::file.path(dir, "cellbrowser.conf")
    base::cat(config, file=confPath)
    cb.dir <- base::paste(dir, "html_data", sep="/")
    cb <- reticulate::import(module = "cellbrowser")
    cb$cellbrowser$build(dir, cb.dir)
}