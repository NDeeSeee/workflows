#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})

suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(future))
suppressMessages(library(DESeq2))
suppressMessages(library(tibble))
suppressMessages(library(scales))
suppressMessages(library(flexmix))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(argparse))
suppressMessages(library(patchwork))
suppressMessages(library(data.table))
suppressMessages(library(reticulate))
suppressMessages(library(sctransform))
suppressMessages(library(rtracklayer))
suppressMessages(library(RColorBrewer))
suppressMessages(library(bestNormalize))
suppressMessages(library(GenomicRanges))
suppressMessages(library(SeuratWrappers))


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


writeSparseMatrix <- function (inMat, outFname, sliceSize=1000){
    fnames <- c()
    mat <- inMat
    geneCount <- nrow(mat)
    startIdx <- 1
    while (startIdx < geneCount){
        endIdx <- min(startIdx + sliceSize - 1, geneCount)
        matSlice <- mat[startIdx:endIdx,]
        denseSlice <- as.matrix(matSlice)
        dt <- data.table(denseSlice)
        dt <- cbind(gene=rownames(matSlice), dt)
        writeHeader <- FALSE
        if (startIdx == 1) { 
            writeHeader <- TRUE
        }
        sliceFname <- paste0("temp", startIdx, ".txt")
        fwrite(dt, sep="\t", file=sliceFname, quote=FALSE, col.names=writeHeader)
        fnames <- append(fnames, sliceFname)
        startIdx <- startIdx + sliceSize
    }
    system(paste("cat", paste(fnames, collapse=" "), "| gzip >", outFname, sep=" "))
    unlink(fnames)
}


findMatrix <- function(object, matrix.slot){
    if (matrix.slot == "counts") {
        counts <- GetAssayData(object=object, slot="counts")
    } else if (matrix.slot == "scale.data") {
        counts <- GetAssayData(object=object, slot="scale.data")
    } else if (matrix.slot == "data") {
        counts <- GetAssayData(object=object)
    } else {
        error("matrix.slot can only be one of: counts, scale.data, data")
    }
}


ExportToCellbrowser <- function(
    object,
    matrix.slot,
    dir,
    cluster.field,
    features,
    meta.fields,
    meta.fields.names
) {
    if (!dir.exists(dir)) {
        dir.create(dir)
    }
    idents <- Idents(object)
    meta <- object@meta.data
    cellOrder <- colnames(object)
    counts <- findMatrix(object, matrix.slot)
    genes <- rownames(x=object)
    dr <- object@reductions

    gzPath <- file.path(dir, "exprMatrix.tsv.gz")
    if ((((ncol(counts)/1000)*(nrow(counts)/1000))>2000) && is(counts, 'sparseMatrix')){
        writeSparseMatrix(counts, gzPath);
    } else {
        mat <- as.matrix(counts)
        df <- as.data.frame(mat, check.names=FALSE)
        df <- data.frame(gene=genes, df, check.names=FALSE)
        z <- gzfile(gzPath, "w")
        write.table(x=df, sep="\t", file=z, quote=FALSE, row.names=FALSE)
        close(con=z)
    }
    embeddings = names(dr)
    embeddings.conf <- c()
    for (embedding in embeddings) {
        emb <- dr[[embedding]]
        df <- emb@cell.embeddings
        if (ncol(df) > 2){
            df <- df[, 1:2]
        }
        colnames(df) <- c("x", "y")
        df <- data.frame(cellId=rownames(df), df, check.names=FALSE)
        write.table(
            df[cellOrder,],
            sep="\t",
            file=file.path(dir, sprintf("%s.coords.tsv", embedding)),
            quote=FALSE,
            row.names=FALSE
        )
        embeddings.conf <- c(
            embeddings.conf,
            sprintf('{"file": "%s.coords.tsv", "shortLabel": "Seurat %1$s"}', embedding)
        )
    }
    df <- data.frame(row.names=cellOrder, check.names=FALSE)
    enum.fields <- c()
    for (i in 1:length(meta.fields)){
        field <- meta.fields[i]
        name <- meta.fields.names[i]
        if (field %in% colnames(meta)){
            df[[name]] <- meta[[field]]
            if (!is.numeric(df[[name]])) {
                enum.fields <- c(enum.fields, name)
            }
        }
    }
    df <- data.frame(Cell=rownames(df), df, check.names=FALSE)
    write.table(
        as.matrix(df[cellOrder,]),
        sep="\t",
        file=file.path(dir, "meta.tsv"),
        quote=FALSE,
        row.names=FALSE
    )
    if (!is.null(features)){
        write.table(
            features,
            sep="\t",
            file=file.path(dir, "quickGenes.tsv"),
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
    enum.string <- paste0("[", paste(paste0('"', enum.fields, '"'), collapse=", "), "]")
    coords.string <- paste0("[", paste(embeddings.conf, collapse = ",\n"), "]")
    config <- sprintf(
        config,
        cluster.field,
        cluster.field,
        enum.string,
        coords.string
    )
    if (!is.null(features)){
        config <- paste(config, 'quickGenesFile="quickGenes.tsv"', sep="\n")
    }
    confPath = file.path(dir, "cellbrowser.conf")
    cat(config, file=confPath)
    cb.dir <- paste(dir, "html_data", sep="/")
    cb <- reticulate::import(module = "cellbrowser")
    cb$cellbrowser$build(dir, cb.dir)
}


export_cellbrowser_data <- function(seurat_data, assay, matrix_slot, resolution, features, rootname){
    tryCatch(
        expr = {
            backup_assay <- DefaultAssay(seurat_data)
            DefaultAssay(seurat_data) <- assay
            meta_fields       <- c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment",       "nucleosome_signal", "frip", "blacklisted_fraction", "miQC.keep")
            meta_fields_names <- c("GEX UMIs",   "Genes",        "Mitochondrial %", "Novelty score",            "ATAC UMIs",   "Peaks",         "TSS enrichment score", "Nucleosome signal", "FRiP", "Blacklisted fraction", "miQC prediction")
            meta_fields <- c(
                meta_fields,
                paste("wsnn_res", resolution, sep=".")
            )
            meta_fields_names <- c(
                meta_fields_names,
                paste("Clustering (", resolution, ")", sep="")
            )

            Idents(seurat_data) <- "condition"
            if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
                meta_fields <- c("condition", meta_fields)
                meta_fields_names <- c("Condition", meta_fields_names)
            }

            cluster_field <- paste("Clustering (", resolution, ")", sep="")[1]
            Idents(seurat_data) <- "new.ident"
            if (length(unique(as.vector(as.character(Idents(seurat_data))))) > 1){
                cluster_field <- "Identity"
                meta_fields <- c("new.ident", meta_fields)
                meta_fields_names <- c("Identity", meta_fields_names)
            }

            tryCatch(
                expr = {
                    custom_fields <- grep("^custom_", colnames(seurat_data@meta.data), value=TRUE, ignore.case=TRUE)
                    custom_fields_names <- gsub("custom_", "", custom_fields)
                    meta_fields <- c(
                        meta_fields,
                        custom_fields
                    )
                    meta_fields_names <- c(
                        meta_fields_names,
                        custom_fields_names
                    )
                },
                error = function(e){
                    print(paste("Failed to add extra metadata fields due to", e))
                }
            )
            ExportToCellbrowser(
                seurat_data,
                matrix.slot=matrix_slot,
                dir=rootname,
                cluster.field=cluster_field,
                features=features,
                meta.fields=meta_fields,
                meta.fields.names=meta_fields_names
            )
            print(paste("Export UCSC Cellbrowser data to", rootname, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export UCSC Cellbrowser data to", rootname, sep=" "))
        },
        finally = {
            DefaultAssay(seurat_data) <- backup_assay
        }
    )
}


extend_metadata_by_barcode <- function(seurat_data, location) {
    tryCatch(
        expr = {
            extra_metadata <- read.table(
                location,
                sep=get_file_type(location),
                header=TRUE,
                check.names=FALSE,
                stringsAsFactors=FALSE
            )
            print(paste("Extra metadata is successfully loaded from ", location))
            refactored_metadata <- data.frame(Cells(seurat_data)) %>%           # create a dataframe with only one column
                                   dplyr::rename("barcode"=1) %>%               # rename that column to barcode
                                   left_join(extra_metadata, by="barcode") %>%  # intersect with loaded extra metadata by "barcode"
                                   remove_rownames() %>%
                                   column_to_rownames("barcode") %>%
                                   replace(is.na(.), "Unknown") %>%             # in case an extra metadata had less barcodes than we had in our Seurat object
                                   rename_with(~paste0("custom_", .x))          # add prefix to all extra metadata columns
            seurat_data <- AddMetaData(
                seurat_data,
                refactored_metadata[Cells(seurat_data), , drop=FALSE]           # to guarantee the proper cells order
            )
        },
        error = function(e){
            print(paste("Failed to add extra cells metadata due to", e))
        },
        finally = {
            return (seurat_data)
        }
    )
}


apply_cell_filters <- function(seurat_data, barcodes_data) {
    filtered_seurat_data <- subset(seurat_data, cells=barcodes_data)
    return (filtered_seurat_data)
}


get_tss_positions <- function(annotation_ranges){
    # Based on GetTSSPositions function from signac/R/utilities.R
    # adapted to work with refgene GTF annotations file
    annotation_df <- as.data.table(x=annotation_ranges)
    annotation_df$strand <- as.character(x=annotation_df$strand)
    annotation_df$strand <- ifelse(
        test = annotation_df$strand == "*",
        yes = "+",
        no = annotation_df$strand
    )
    collapsed_annotation_df <- annotation_df[
        , .(unique(seqnames),
            min(start),
            max(end),
            strand[[1]],
            gene_name[[1]]),
        "gene_id"
    ]
    colnames(x=collapsed_annotation_df) <- c(
        "gene_id", "seqnames", "start", "end", "strand", "gene_name"
    )
    collapsed_annotation_df$gene_name <- make.unique(names=collapsed_annotation_df$gene_name)
    collapsed_ranges <- makeGRangesFromDataFrame(
        df=collapsed_annotation_df,
        keep.extra.columns=TRUE
    )
    tss_positions <- resize(collapsed_ranges, width=1, fix="start")
    return (tss_positions)
}


add_qc_metrics <- function(seurat_data, blacklisted_data, args) {
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "RNA"
    seurat_data$log10_gene_per_log10_umi <- log10(seurat_data$nFeature_RNA) / log10(seurat_data$nCount_RNA)
    seurat_data$mito_percentage <- PercentageFeatureSet(seurat_data, pattern=args$mitopattern) 
    if (args$skipmiqc){
        print(
            paste(
                "Skipping threshold prediction for the percentage of",
                "transcripts mapped to mitochondrial genes (do not run MiQC)"
            )
        )
    } else {
        print(
            paste(
                "Trying to predict the threshold for the percentage of transcripts",
                "mapped to mitochondrial genes using MiQC"
            )
        )
        tryCatch(
            expr = {
                Idents(seurat_data) <- "new.ident"                                      # safety measure
                identities <- unique(as.vector(as.character(Idents(seurat_data))))
                merged_seurat_data <- NULL
                for (i in 1:length(identities)) {
                    current_identity <- identities[i]
                    filtered_seurat_data <- subset(seurat_data, idents=current_identity)
                    filtered_seurat_data <- RunMiQC(
                        filtered_seurat_data,
                        percent.mt="mito_percentage",
                        nFeature_RNA="nFeature_RNA",
                        model.type="linear",
                        posterior.cutoff=0.75,
                        model.slot="flexmix_model",
                        verbose=args$verbose
                    )
                    if (is.null(merged_seurat_data)){
                        merged_seurat_data <- filtered_seurat_data
                    } else {
                        merged_seurat_data <- merge(merged_seurat_data, y=filtered_seurat_data)
                    }
                }
                seurat_data <- merged_seurat_data
            },
            error = function(e){
                print(paste(" ", "Failed running MiQC due to", e))
            }
        )
    }

    DefaultAssay(seurat_data) <- "ATAC"

    # Nucleosome banding pattern
    #   The histogram of DNA fragment sizes (determined from the paired-end sequencing reads)
    #   should exhibit a strong nucleosome banding pattern corresponding to the length of DNA
    #   wrapped around a single nucleosome. We calculate this per single cell, and quantify
    #   the approximate ratio of mononucleosomal to nucleosome-free fragments (stored as 
    #   nucleosome_signal)
    seurat_data <- NucleosomeSignal(seurat_data, verbose=args$verbose)

    # Transcriptional start site (TSS) enrichment score
    #   The ENCODE project has defined an ATAC-seq targeting score based on the ratio of fragments
    #   centered at the TSS to fragments in TSS-flanking regions (see https://www.encodeproject.org/data-standards/terms/).
    #   Poor ATAC-seq experiments typically will have a low TSS enrichment score. We can compute this
    #   metric for each cell with the TSSEnrichment() function, and the results are stored in metadata
    #   under the column name TSS.enrichment.
    #   We need to calculate tss.positions by ourselves as those that are calculated automatically
    #   will require biotypes="protein_coding" filed in GTF file
    #   see https://m.ensembl.org/info/genome/genebuild/biotypes.html for details about biotypes in GTF
    seurat_data <- TSSEnrichment(
        seurat_data,
        tss.positions=get_tss_positions(Annotation(seurat_data[["ATAC"]])),
        fast=FALSE,                                                          # set fast=FALSE, because we want to build TSS Enrichment plot later
        verbose=args$verbose
    )

    # Fraction of fragments in peaks.
    #   Represents the fraction of all fragments that fall within ATAC-seq peaks. Cells with low values
    #   (i.e. <15-20%) often represent low-quality cells or technical artifacts that should be removed.
    #   Note that this value can be sensitive to the set of peaks used.
    fragments_data <- CountFragments(
        fragments=args$fragments,
        cells=colnames(seurat_data),           # limit it to only those cells that are present in seurat_data
        verbose=args$verbose
    ) %>% column_to_rownames("CB")             # for easy access to cell barcodes
    seurat_data$fragments <- fragments_data[colnames(seurat_data), "frequency_count"]  # select by rownames to make sure the cells order wasn't accidentally changed
    seurat_data <- FRiP(
        seurat_data,
        assay="ATAC",                          # FRiP can't take the default assay, so we set it explicitly
        total.fragments="fragments",
        col.name="frip",
        verbose=args$verbose
    )

    # Ratio reads in genomic blacklist regions.
    #   The ENCODE project has provided a list of blacklist regions, representing reads which are often
    #   associated with artefactual signal. Cells with a high proportion of reads mapping to these areas
    #   often represent technical artifacts and should be removed.
    if (!is.null(blacklisted_data)){
        seurat_data$blacklisted_fraction <- FractionCountsInRegion(
            seurat_data, 
            assay="ATAC",
            regions=blacklisted_data
        )
    } else {
        seurat_data$blacklisted_fraction <- 0  # blacklisted regions file wasn't provided, so we set everything to 0
    }
    DefaultAssay(seurat_data) <- backup_assay
    return (seurat_data)
}


integrate_atac_data <- function(seurat_data, args) {
    DefaultAssay(seurat_data) <- "ATAC"                                          # safety measure in case ATAC wasn't the default assay
    splitted_seurat_data <- SplitObject(seurat_data, split.by="new.ident")
    if (args$skipatacntrg | length(splitted_seurat_data) == 1){
        print("Skipping datasets integration: either forced or only one identity is present")
        scaled_norm_seurat_data <- FindTopFeatures(splitted_seurat_data[[1]], min.cutoff=args$highvaratac, verbose=args$verbose)
        scaled_norm_seurat_data <- RunTFIDF(scaled_norm_seurat_data, verbose=args$verbose)
        scaled_norm_seurat_data <- RunSVD(                                                          # by default computes 50 singular values
            scaled_norm_seurat_data,
            verbose=args$verbose,
            reduction.name="atac_lsi"
        )
        return (scaled_norm_seurat_data)
    } else {
        print("Running ATAC datasets integration using IntegrateEmbeddings algorithm")
        backup_pca_reduction <- seurat_data[["pca"]]          # need to backup reductions generated from RNA assay as the will be deleted by IntegrateEmbeddings function
        backup_rnaumap_reduction <- seurat_data[["rnaumap"]]
        seurat_data <- FindTopFeatures(seurat_data, min.cutoff=args$highvaratac, verbose=args$verbose)
        seurat_data <- RunTFIDF(seurat_data, verbose=args$verbose)
        seurat_data <- RunSVD(seurat_data, verbose=args$verbose)                                    # by default computes 50 singular values
        for (i in 1:length(splitted_seurat_data)) {
            splitted_seurat_data[[i]] <- FindTopFeatures(splitted_seurat_data[[i]], min.cutoff=args$highvaratac, verbose=args$verbose)
            splitted_seurat_data[[i]] <- RunTFIDF(splitted_seurat_data[[i]], verbose=args$verbose)
            splitted_seurat_data[[i]] <- RunSVD(splitted_seurat_data[[i]], verbose=args$verbose)    # by default computes 50 singular values
        }
        integration_anchors <- FindIntegrationAnchors(
            splitted_seurat_data,
            anchor.features=rownames(seurat_data),
            reduction="rlsi",
            dims=2:50,                                                                              # use up to computed singular values
            verbose=args$verbose
        )
        integrated_seurat_data <- IntegrateEmbeddings(
            anchorset=integration_anchors,
            reductions=seurat_data[["lsi"]],
            new.reduction.name="atac_lsi",
            k.weight=min(min(table(Idents(seurat_data))), 100),  # 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering
            dims.to.integrate=1:50                                                                  # use up to computed singular values
        )
        integrated_seurat_data[["pca"]] <- backup_pca_reduction                                     # restoring deleted reductions
        integrated_seurat_data[["rnaumap"]] <- backup_rnaumap_reduction
        return (integrated_seurat_data)
    }
}


integrate_gex_data <- function(seurat_data, args) {
    DefaultAssay(seurat_data) <- "RNA"                                                       # safety measure, in case RNA was not default assay
    splitted_seurat_data <- SplitObject(seurat_data, split.by="new.ident")
    vars_to_regress <- NULL                                                                  # used in all three places, so no need to repeat it each time
    if (args$regressmt) {
        vars_to_regress <- c("mito_percentage")
    }
    if (args$skipgexntrg | length(splitted_seurat_data) == 1){
        print(
            paste(
                "Skipping datasets integration: either forced or only one identity is present.",
                "Running log-normalization and scalling instead."
            )
        )
        scaled_norm_seurat_data <- NormalizeData(splitted_seurat_data[[1]], verbose=args$verbose)
        scaled_norm_seurat_data <- FindVariableFeatures(
            scaled_norm_seurat_data,
            nfeatures=args$highvargex,
            verbose=args$verbose
        )
        scaled_norm_seurat_data <- ScaleData(
            scaled_norm_seurat_data,
            vars.to.regress=vars_to_regress,
            verbose=args$verbose
        )
        DefaultAssay(scaled_norm_seurat_data) <- "RNA"
        return (scaled_norm_seurat_data)
    } else if (args$nosct) {
        print("Skipping SCTransform and running standard integration algorithm")
        for (i in 1:length(splitted_seurat_data)) {
            splitted_seurat_data[[i]] <- NormalizeData(splitted_seurat_data[[i]], verbose=args$verbose)
            splitted_seurat_data[[i]] <- FindVariableFeatures(
                splitted_seurat_data[[i]],
                nfeatures=args$highvargex,
                verbose=args$verbose
            )
        }
        integration_features <- SelectIntegrationFeatures(
            splitted_seurat_data,
            nfeatures=args$highvargex,
            verbose=args$verbose
        )
        integration_anchors <- FindIntegrationAnchors(
            splitted_seurat_data,
            normalization.method="LogNormalize",
            anchor.features=integration_features,
            verbose=args$verbose
        )
        integrated_seurat_data <- IntegrateData(
            integration_anchors, 
            new.assay.name="gex_integrated",
            normalization.method="LogNormalize",
            k.weight=min(min(table(Idents(seurat_data))), 100),  # 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering
            verbose=args$verbose
        )
        integrated_seurat_data <- ScaleData(
            integrated_seurat_data,
            vars.to.regress=vars_to_regress,
            verbose=args$verbose
        )
        DefaultAssay(integrated_seurat_data) <- "gex_integrated"
        return (integrated_seurat_data)
    } else {
        print("Running SCTransform integration algorithm")
        for (i in 1:length(splitted_seurat_data)) {
            splitted_seurat_data[[i]] <- SCTransform(
                splitted_seurat_data[[i]],
                assay="RNA",
                new.assay.name="SCT",
                variable.features.n=args$highvargex,
                vars.to.regress=vars_to_regress,
                verbose=args$verbose
            )
        }
        integration_features <- SelectIntegrationFeatures(
            splitted_seurat_data,
            nfeatures=args$highvargex,
            verbose=args$verbose
        )
        splitted_seurat_data <- PrepSCTIntegration(
            splitted_seurat_data, 
            anchor.features=integration_features,
            verbose=args$verbose
        )
        integration_anchors <- FindIntegrationAnchors(
            splitted_seurat_data,
            normalization.method="SCT",
            anchor.features=integration_features,
            verbose=args$verbose
        )
        integrated_seurat_data <- IntegrateData(
            integration_anchors, 
            new.assay.name="gex_integrated",
            normalization.method="SCT",
            k.weight=min(min(table(Idents(seurat_data))), 100),  # 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering
            verbose=args$verbose
        )
        DefaultAssay(integrated_seurat_data) <- "gex_integrated"
        return (integrated_seurat_data)
    }
}


get_all_putative_markers <- function(seurat_data, args, assay="RNA", min_diff_pct=-Inf){
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    all_putative_markers <- NULL
    for (i in 1:length(args$resolution)) {
        resolution <- args$resolution[i]
        Idents(seurat_data) <- paste("wsnn_res", resolution, sep=".")
        tryCatch(
            expr = {
                markers <- FindAllMarkers(
                    seurat_data,
                    logfc.threshold=args$gexlogfc,
                    min.pct=args$gexminpct,
                    only.pos=args$gexonlypos,
                    test.use=args$gextestuse,
                    min.diff.pct=min_diff_pct,
                    verbose=FALSE
                ) %>% relocate(cluster, gene, .before=1)
                if (nrow(markers) > 0) {
                    markers <- markers %>% cbind(resolution=resolution, .)
                } else {
                    markers <- markers %>% add_column(resolution=numeric(), .before=1)  # safety measure in case markers was empty
                }
                if (!is.null(all_putative_markers)) {
                    all_putative_markers <- rbind(all_putative_markers, markers)
                } else {
                    all_putative_markers <- markers
                }
            },
            error = function(e){
                print(paste("Failed find putative gene markers for resolution", resolution, "due to", e))
            }
        )
        Idents(seurat_data) <- "new.ident"
    }
    DefaultAssay(seurat_data) <- backup_assay
    return (all_putative_markers)
}


apply_qc_filters <- function(seurat_data, cell_identity_data, args) {
    merged_seurat_data <- NULL
    for (i in 1:length(args$mingenes)){
        identity <- cell_identity_data$library_id[i]

        mingenes <- args$mingenes[i]
        maxgenes <- args$maxgenes[i]
        gexminumi <- args$gexminumi[i]
        minnovelty <- args$minnovelty[i]
        atacminumi <- args$atacminumi[i]
        maxnuclsignal <- args$maxnuclsignal[i]
        mintssenrich <- args$mintssenrich[i]
        minfrip <- args$minfrip[i]
        maxblacklisted <- args$maxblacklisted[i]

        print(paste("Filtering", identity))
        print(paste(" ", mingenes, "<= Genes per cell <=", maxgenes))
        print(paste(" ", "GEX UMIs per cell >=", gexminumi))
        print(paste(" ", "GEX novelty score >=", minnovelty))
        print(paste(" ", "ATAC UMIs per cell >=", atacminumi))
        print(paste(" ", "Nucleosome signal <=", maxnuclsignal))
        print(paste(" ", "TSS enrichment score >=", mintssenrich))
        print(paste(" ", "FRiP >=", minfrip))
        print(paste(" ", "Ratio of reads in genomic blacklist regions <=", maxblacklisted))
        print(paste(" ", "Percentage of GEX transcripts mapped to mitochondrial genes <=", args$maxmt))

        filtered_seurat_data <- subset(
            seurat_data,
            idents=identity,
            subset=(nFeature_RNA >= mingenes) &
                (nFeature_RNA <= maxgenes) &
                (nCount_RNA >= gexminumi) &
                (log10_gene_per_log10_umi >= minnovelty) &
                (nCount_ATAC >= atacminumi) &
                (nucleosome_signal <= maxnuclsignal) &
                (TSS.enrichment >= mintssenrich) &
                (frip >= minfrip) &
                (blacklisted_fraction <= maxblacklisted) &
                (mito_percentage <= args$maxmt)
        )

        if (is.null(merged_seurat_data)){
            merged_seurat_data <- filtered_seurat_data
        } else {
            merged_seurat_data <- merge(merged_seurat_data, y=filtered_seurat_data)
        }
    }
    return (merged_seurat_data)
}


export_geom_bar_plot <- function(data, rootname, x_axis, color_by, x_label, y_label, legend_title, plot_title, palette="Paired", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                geom_bar(colour="black") +
                geom_text(stat="count", aes(label=..count..), vjust=-1) +
                xlab(x_label) +
                ylab(y_label) +
                guides(fill=guide_legend(legend_title), x=guide_axis(angle=45)) +
                ggtitle(plot_title) +
                scale_fill_brewer(palette=palette)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export geom bar plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export geom bar plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_geom_density_plot <- function(data, rootname, x_axis, color_by, facet_by, x_left_intercept, x_label, y_label, legend_title, plot_title, x_right_intercept=NULL,  scale_x_log10=FALSE, scale_y_log10=FALSE, zoom_on_intercept=FALSE, alpha=0.9, show_ranked=FALSE, ranked_x_label="Ranked cells", palette="Paired", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            palette_colors <- RColorBrewer::brewer.pal(99, palette)  # use 99 to get all available colors. As we added LETTERS, the order will be correct
            intercept_data <- data %>%
                              dplyr::select(all_of(color_by), all_of(facet_by)) %>%
                              distinct() %>%
                              arrange(all_of(color_by)) %>%
                              add_column(color=palette_colors[1:nrow(.)],x_left=x_left_intercept)

            plot <- ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                    geom_density(alpha=alpha) +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(fill=guide_legend(legend_title)) +
                    ggtitle(plot_title) +
                    facet_wrap(as.formula(paste("~", facet_by))) +
                    scale_fill_brewer(palette=palette) +
                    geom_vline(intercept_data, mapping=aes(xintercept=x_left), color=intercept_data$color, alpha=0.7) +
                    geom_label_repel(
                        intercept_data, mapping=aes(x=x_left, y=Inf, label=x_left),
                        color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                        show.legend=FALSE
                    )

            if (!is.null(x_right_intercept)){
                intercept_data <- intercept_data %>% add_column(x_right=x_right_intercept)
                plot <- plot +
                        geom_vline(intercept_data, mapping=aes(xintercept=x_right), color=intercept_data$color, alpha=0.7) +
                        geom_label_repel(
                            intercept_data, mapping=aes(x=x_right, y=Inf, label=x_right),
                            color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                            show.legend=FALSE
                        )
            }

            if (scale_x_log10){ plot <- plot + scale_x_log10() }
            if (scale_y_log10){ plot <- plot + scale_y_log10() }

            if (zoom_on_intercept) {
                zoomed_plot <- ggplot(data, aes_string(x=x_axis, color=color_by)) +
                               geom_density(show.legend=FALSE, size=2) +
                               xlab(x_label) +
                               ylab(y_label) +
                               scale_color_brewer(palette=palette) +
                               geom_vline(intercept_data, mapping=aes(xintercept=x_left), color=intercept_data$color, alpha=0.7) +
                               geom_label_repel(
                                   intercept_data, mapping=aes(x=x_left, y=Inf, label=x_left),
                                   color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                                   show.legend=FALSE
                               )
                if (!is.null(x_right_intercept)){
                    zoomed_plot <- zoomed_plot +
                                   coord_cartesian(xlim=c(min(intercept_data$x_left), max(intercept_data$x_right))) +
                                   geom_vline(intercept_data, mapping=aes(xintercept=x_right), color=intercept_data$color, alpha=0.7) +
                                   geom_label_repel(
                                       intercept_data, mapping=aes(x=x_right, y=Inf, label=x_right),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="y", size=3,
                                       show.legend=FALSE
                                   )
                } else {
                    zoomed_plot <- zoomed_plot + coord_cartesian(xlim=c(NA, max(intercept_data$x_left)))
                }
                if (scale_x_log10){ zoomed_plot <- zoomed_plot + scale_x_log10() }
                if (scale_y_log10){ zoomed_plot <- zoomed_plot + scale_y_log10() }
                plot <- plot / zoomed_plot
            }

            if (show_ranked) {
                ranked_plot <- data %>%
                               arrange(get(color_by), get(x_axis)) %>%
                               ggplot(aes_string(x=seq_along(data[[x_axis]]), y=x_axis, color=color_by)) +
                               geom_point(show.legend=FALSE, size=0.5) +
                               xlab(ranked_x_label) +
                               ylab(x_label) +
                               scale_y_log10() +
                               scale_color_brewer(palette=palette) +
                               geom_hline(intercept_data, mapping=aes(yintercept=x_left), color=intercept_data$color, alpha=0.7) +
                               geom_label_repel(
                                   intercept_data, mapping=aes(x=Inf, y=x_left, label=x_left),
                                   color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="x", size=3,
                                   show.legend=FALSE
                               )
                if (!is.null(x_right_intercept)){
                    ranked_plot <- ranked_plot +
                                   geom_hline(intercept_data, mapping=aes(yintercept=x_right), color=intercept_data$color, alpha=0.7) +
                                   geom_label_repel(
                                       intercept_data, mapping=aes(x=-Inf, y=x_right, label=x_right),
                                       color="black", fill=intercept_data$color, alpha=0.7, segment.colour=NA, direction="x", size=3,
                                       show.legend=FALSE
                                   )
                }
                plot <- plot / ranked_plot
            }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export geom density plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export geom density plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_geom_point_plot <- function(data, rootname, x_axis, y_axis, facet_by, x_left_intercept, y_low_intercept, color_by, colors, color_limits, color_break, x_label, y_label, legend_title, plot_title, y_high_intercept=NULL, scale_x_log10=FALSE, scale_y_log10=FALSE, alpha=0.2, alpha_intercept=0.5, palette="Paired", pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            palette_colors <- RColorBrewer::brewer.pal(99, palette)  # use 99 to get all available colors. As we added LETTERS, the order will be correct
            intercept_data <- data %>%
                              dplyr::select(all_of(facet_by)) %>%
                              distinct() %>%
                              arrange(all_of(facet_by)) %>%
                              add_column(color=palette_colors[1:nrow(.)], x_left=x_left_intercept, y_low=y_low_intercept)

            if (!is.null(y_high_intercept)){
                intercept_data <- intercept_data %>% add_column(y_high=y_high_intercept)
            }

            plot <- ggplot(data, aes_string(x=x_axis, y=y_axis, color=color_by)) +
                    geom_point(alpha=alpha) +
                    scale_colour_gradientn(
                        colours=c(colors[1], colors),
                        values=rescale(c(color_limits[1], color_break-0.01*color_break, color_break, color_limits[2])),
                        breaks=c(color_break),
                        limits=color_limits
                    ) +
                    stat_smooth(method=lm) +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(color=guide_colourbar(legend_title)) +
                    ggtitle(plot_title) +
                    facet_wrap(as.formula(paste("~", facet_by))) +
                    geom_vline(intercept_data, mapping=aes(xintercept=x_left), color=intercept_data$color, alpha=alpha_intercept) +
                    geom_hline(intercept_data, mapping=aes(yintercept=y_low), color=intercept_data$color, alpha=alpha_intercept) +
                    geom_label_repel(
                        intercept_data, mapping=aes(x=x_left, y=Inf, label=x_left),
                        color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="y", size=3,
                        show.legend=FALSE
                    ) +
                    geom_label_repel(
                        intercept_data, mapping=aes(x=Inf, y=y_low, label=y_low),
                        color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="x", size=3,
                        show.legend=FALSE
                    )

            if (!is.null(y_high_intercept)){
                plot <- plot +
                        geom_hline(intercept_data, mapping=aes(yintercept=y_high), color=intercept_data$color, alpha=alpha_intercept) +
                        geom_label_repel(
                            intercept_data, mapping=aes(x=Inf, y=y_high, label=y_high),
                            color="black", fill=intercept_data$color, alpha=alpha_intercept, direction="x", size=3,
                            show.legend=FALSE
                        )
            }

            if (scale_x_log10){ plot <- plot + scale_x_log10() }
            if (scale_y_log10){ plot <- plot + scale_y_log10() }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export geom point plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export geom point plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_feature_scatter_plot <- function(data, rootname, x_axis, y_axis, x_label, y_label, split_by, color_by, plot_title, legend_title, combine_guides=NULL, palette=NULL, alpha=NULL, jitter=FALSE, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            Idents(data) <- split_by
            identities <- unique(as.vector(as.character(Idents(data))))
            plots <- list()
            for (i in 1:length(identities)) {
                current_identity <- identities[i]
                filtered_data <- subset(data, idents=current_identity)
                plots[[current_identity]] <- FeatureScatter(
                    filtered_data,
                    feature1=x_axis,
                    feature2=y_axis,
                    group.by=color_by,
                    plot.cor=FALSE,         # will be overwritten by title anyway
                    jitter=jitter
                )
            }
            Idents(data) <- "new.ident"
            plots <- lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] +
                              ggtitle(identities[i]) +
                              xlab(x_label) +
                              ylab(y_label) +
                              guides(color=guide_legend(legend_title)) +
                              theme_gray()
                if (!is.null(palette)) { plots[[i]] <- plots[[i]] + scale_color_manual(values=palette) }
                if (!is.null(alpha)) { plots[[i]]$layers[[1]]$aes_params$alpha <- alpha }
                return (plots[[i]])
            })
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export feature scatter plot to ", rootname, ".(png/pdf)", sep=""))

        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export feature scatter plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
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
            print(paste("Failed to export data to", location, sep=" "))
        }
    )
}


get_density_data <- function(data, features, labels, assay="RNA", slot="data", n_bins=1024, x_min=NULL, x_max=NULL){
    DefaultAssay(data) <- assay                               # no need to backup assay as we don't return data from the function
    slot_data <- FetchData(data, vars=features, slot=slot)
    density_data <- as.data.frame(
        do.call(
            cbind,                                            # merges to one big dataframe
            lapply(                                           # returns list of dataframes
                features,
                function(f){
                    feature_data <- slot_data[, f]
                    l <- labels[match(f, features)]
                    density_ls <- density(
                        feature_data,
                        from=ifelse(is.null(x_min), 0, x_min),
                        to=ifelse(is.null(x_max), max(feature_data[is.finite(feature_data)]), x_max),
                        n=n_bins
                    )
                    density_df <- data.frame(density_ls$x, density_ls$y)
                    colnames(density_df) <- c(paste(l, "x"), paste(l, "y"))
                    return (density_df)
                }
            )
        )
    )
    return (density_data)
}


export_vln_plot <- function(data, features, labels, rootname, plot_title, legend_title, log=FALSE, group_by=NULL, hide_x_text=FALSE, pt_size=NULL, palette=NULL, combine_guides=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plots <- VlnPlot(
                         data,
                         features=features,
                         pt.size=pt_size,
                         group.by=group_by,
                         log=log,
                         combine=FALSE       # to return a list of gglots
                     )
            plots <- lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] +
                              ggtitle(labels[i]) +
                              theme_gray() +
                              theme(axis.title.x=element_blank()) +
                              guides(fill=guide_legend(legend_title)) +
                              stat_boxplot(width=0.15, geom="errorbar") +
                              geom_boxplot(width=0.15, outlier.alpha=0) +
                              RotatedAxis()
                if (!is.null(palette)){ plots[[i]] <- plots[[i]] + scale_fill_brewer(palette=palette) }
                if (hide_x_text){ plots[[i]] <- plots[[i]] + theme(axis.text.x=element_blank()) }
                return (plots[[i]])
            })
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export violin plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export violin plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_dim_plot <- function(data, rootname, reduction, plot_title, legend_title, split_by=NULL, group_by=NULL, perc_split_by=NULL, perc_group_by=NULL, label=FALSE, palette=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- DimPlot(
                        data,
                        reduction=reduction,
                        split.by=split_by,
                        group.by=group_by,
                        label=label
                    ) +
                    theme_gray() +
                    ggtitle(plot_title) +
                    guides(color=guide_legend(legend_title, override.aes=list(size=3)))

            if (!is.null(palette)){ plot <- plot + scale_color_brewer(palette=palette) }

            if(!is.null(perc_split_by) && !is.null(perc_group_by)){
                width <- 2 * width
                perc_data <- data@meta.data %>%
                             dplyr::group_by(across(all_of(perc_split_by)), across(all_of(perc_group_by))) %>%
                             summarise(counts=n(), .groups="drop_last") %>%               # drop the perc_group_by grouping level so we can get only groups defined by perc_split_by
                             mutate(freq=counts/sum(counts)*100) %>%                      # sum is taken for the group defined by perc_split_by
                             arrange(all_of(perc_split_by), all_of(perc_group_by))        # sort for consistency
                label_data <- data@meta.data %>%
                              dplyr::group_by(across(all_of(perc_split_by))) %>%
                              summarise(counts=n(), .groups="drop_last") %>%              # drops all grouping as we have only one level
                              arrange(all_of(perc_split_by))                              # sort for consistency
                perc_plot <- ggplot(perc_data, aes_string(x=perc_split_by, y="freq", fill=perc_group_by)) +
                             geom_col(position="dodge", width=0.9, linetype="solid", color="black", show.legend=FALSE) +
                             xlab("") +
                             ylab("Cells percentage") +
                             theme_gray() +
                             geom_label_repel(
                                label_data, mapping=aes(y=-Inf, label=counts),
                                color="black", fill="white", segment.colour=NA,
                                direction="y", size=3, show.legend=FALSE
                             ) +
                             RotatedAxis()
                plot <- plot + perc_plot
            }

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export dim plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export dim plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_elbow_plot <- function(data, rootname, plot_title, ndims=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- ElbowPlot(
                        data,
                        ndims=ndims,
                        reduction="pca"
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

            print(paste("Export elbow plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export elbow plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_corr_plot <- function(data, reduction, qc_columns, qc_labels, plot_title, rootname, ndims=NULL, combine_guides=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            embeddings <- Embeddings(data[[reduction]])
            ndims=ifelse(is.null(ndims), length(data[[reduction]]), ndims)
            plots <- list()
            for (i in 1:length(qc_columns)) {
                current_qc_column <- qc_columns[i]
                current_qc_label <- qc_labels[i]
                if ( !(qc_columns[i] %in% colnames(data@meta.data)) ){
                    print(
                        paste(
                            "Column", current_qc_column, "was not found,",
                            "skipping", current_qc_label
                        )
                    )
                    next
                }
                qc_data <- data[[current_qc_column]]
                corr_data <- as.data.frame(cor(x=embeddings, y=qc_data))
                corr_data$correlation <- corr_data[, 1]
                corr_data$dimension <- seq_len(length.out=nrow(corr_data))

                plots[[current_qc_column]] <- ggplot(corr_data, aes(dimension, correlation)) +
                                              geom_point() +
                                              xlab("Dimension") +
                                              ylab("Correlation") +
                                              xlim(c(0, ndims)) +
                                              ylim(c(-1, 1)) +
                                              theme_gray() +
                                              ggtitle(current_qc_label)
            }
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export correlation plot to ", rootname, ".(png/pdf)", sep=""))

        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export correlation plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_tss_plot <- function(data, rootname, plot_title, split_by, group_by_value=NULL, combine_guides=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            Idents(data) <- split_by
            identities <- unique(as.vector(as.character(Idents(data))))
            plots <- list()
            for (i in 1:length(identities)) {
                current_identity <- identities[i]
                filtered_data <- subset(data, idents=current_identity)
                group_by <- NULL
                if (!is.null(group_by_value)){
                    filtered_data$tss_group_by <- ifelse(
                        filtered_data$TSS.enrichment >= group_by_value ,
                        paste("High ", "(bigger or equal to ", group_by_value, ")", sep=""),
                        paste("Low ", "(smaller than ", group_by_value, ")", sep="")
                    )
                    group_by <- "tss_group_by"
                }
                plots[[current_identity]] <- TSSPlot(
                        filtered_data,
                        group.by=group_by
                    ) +
                    ggtitle(current_identity) +
                    theme_gray() +
                    NoLegend()
            }
            Idents(data) <- "new.ident"
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export TSS Enrichment plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export TSS Enrichment plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_fragments_hist <- function(data, rootname, plot_title, split_by, group_by_value=NULL, combine_guides=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            Idents(data) <- split_by
            identities <- unique(as.vector(as.character(Idents(data))))
            plots <- list()
            for (i in 1:length(identities)) {
                current_identity <- identities[i]
                filtered_data <- subset(data, idents=current_identity)
                group_by <- NULL
                if (!is.null(group_by_value)){
                    filtered_data$ns_group_by <- ifelse(
                        filtered_data$nucleosome_signal >= group_by_value ,
                        paste("High nucleosome signal ", "(bigger or equal to ", group_by_value, ")", sep=""),
                        paste("Low nucleosome signal", "(smaller than ", group_by_value, ")", sep="")
                    )
                    group_by <- "ns_group_by"
                }
                plots[[current_identity]] <- FragmentHistogram(
                        filtered_data,
                        group.by=group_by
                    ) +
                    theme_gray() +
                    ggtitle(current_identity) +
                    NoLegend()
            }
            Idents(data) <- "new.ident"
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export fragments length histogram to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called dev.off() with error -", e))})
            print(paste("Failed to export fragments length histogram to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


get_qc_metrics_pca <- function(seurat_data, qc_columns, qc_labels, orq_transform=FALSE){
    tryCatch(
        expr = {
            print("Computing PCA for the following QC metrics")
            print(paste(qc_labels, collapse=", "))
            qc_columns_corrected <- c()
            qc_labels_corrected <- c()
            for (i in 1:length(qc_columns)){
                if (qc_columns[i] %in% colnames(seurat_data@meta.data)){
                    qc_columns_corrected <- c(qc_columns_corrected, qc_columns[i])
                    qc_labels_corrected <- c(qc_labels_corrected, qc_labels[i])
                } else {
                    print(
                        paste(
                            "Column", qc_columns[i], "was not found,",
                            "skipping", qc_labels[i]
                        )
                    )
                }
            }
            target_data <- as.data.frame(seurat_data[[qc_columns_corrected]]) %>% 
                           filter_all(all_vars(!is.infinite(.)))                          # safety measure
            print(
                paste(
                    "Rows removed due to having infinite values",
                    "in any of the selected column -",
                    nrow(seurat_data@meta.data) - nrow(target_data)
                )
            )
            if (orq_transform){
                print("Running Ordered Quantile (ORQ) normalization transformation")
                target_data <- target_data %>% mutate_all(function(x){return (orderNorm(x)$x.t)})
            }
            pca_raw <- prcomp(
                t(target_data),
                center=!orq_transform,         # no need to center or scale when data is already ORQ-transformed
                scale.=!orq_transform
            )
            print(summary(pca_raw))
            pca_scores <- as.data.frame(pca_raw$x)
            pca_scores$labels <- qc_labels_corrected
            pca_variance <- round(pca_raw$sdev / sum(pca_raw$sdev) * 100, 2)
            return (list(scores=pca_scores, variance=pca_variance))
        },
        error = function(e){
            print(paste("Failed to compute PCA for QC metrics due to", e))
        }
    )
}


export_pca_plot <- function(pca_data, pcs, rootname, plot_title, legend_title, color_by="label", label_size=5, pt_size=8, pt_shape=19, alpha=0.75, palette="Paired", pdf=FALSE, width=1200, height=1200, resolution=100){
    tryCatch(
        expr = {
            x_score_column <- paste0("PC", pcs[1])
            y_score_column <- paste0("PC", pcs[2])
            x_variance <- pca_data$variance[pcs[1]]
            y_variance <- pca_data$variance[pcs[2]]
            plot <- ggplot(
                        pca_data$scores,
                        aes_string(x=x_score_column, y=y_score_column, color=color_by)
                    ) +
                    geom_point(size=pt_size, shape=pt_shape, alpha=alpha) +
                    xlab(paste0(x_score_column, ": ", x_variance, "% variance")) +
                    ylab(paste0(y_score_column, ": ", y_variance, "% variance")) + 
                    geom_label_repel(
                        aes_string(label=color_by),
                        size=label_size,
                        point.padding=0.5,
                        box.padding=0.5,
                        check_overlap=TRUE,
                        show.legend=FALSE
                    ) +
                    ggtitle(plot_title) +
                    guides(color=guide_legend(legend_title)) +
                    scale_color_brewer(palette=palette) +
                    theme_gray()

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export PCA plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export PCA plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_all_qc_plots <- function(seurat_data, suffix, args){
    qc_metrics_pca <- get_qc_metrics_pca(
        seurat_data=seurat_data,
        qc_columns=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal", "frip", "blacklisted_fraction"),
        qc_labels=c("GEX UMIs", "Genes", "Mitochondrial %", "Novelty score", "ATAC UMIs", "Peaks", "TSS enrichment score", "Nucleosome signal", "FRiP", "Fr. of reads in bl-ted reg."),
        orq_transform=TRUE
    )
    export_pca_plot(
        pca_data=qc_metrics_pca,
        rootname=paste(args$output, suffix, "pca", paste(c(1, 2) ,collapse="_"), "qc_mtrcs", sep="_"),
        pcs=c(1, 2),
        plot_title=paste(
            paste(
                paste0("PC", c(1, 2)),
                collapse=" and "
            ),
            " of ORQ-transformed QC metrics PCA (", suffix, ")", sep=""
        ),
        legend_title="QC metrics",
        color_by="labels",
        pdf=args$pdf
    )
    export_pca_plot(
        pca_data=qc_metrics_pca,
        rootname=paste(args$output, suffix, "pca", paste(c(2, 3) ,collapse="_"), "qc_mtrcs", sep="_"),
        pcs=c(2, 3),
        plot_title=paste(
            paste(
                paste0("PC", c(2, 3)),
                collapse=" and "
            ),
            " of ORQ-transformed QC metrics PCA (", suffix, ")", sep=""
        ),
        legend_title="QC metrics",
        color_by="labels",
        pdf=args$pdf
    )
    export_geom_bar_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "cell_count", sep="_"),
        x_axis="new.ident",
        color_by="new.ident",
        x_label="Identity",
        y_label="Cells",
        legend_title="Identity",
        plot_title=paste("Number of cells per dataset (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gex_umi_dnst", sep="_"),
        x_axis="nCount_RNA",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$gexminumi,
        x_label="UMIs per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("GEX UMI density per cell (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "atac_umi_dnst", sep="_"),
        x_axis="nCount_ATAC",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$atacminumi,
        x_label="UMIs per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("ATAC UMI density per cell (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gene_dnst", sep="_"),
        x_axis="nFeature_RNA",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$mingenes,
        x_right_intercept=args$maxgenes,
        x_label="Genes per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Gene density per cell (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "peak_dnst", sep="_"),
        x_axis="nFeature_ATAC",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=0,
        # x_right_intercept=args$maxgenes,
        x_label="Peaks per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Peak density per cell (", suffix, ")", sep=""),
        scale_x_log10=FALSE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "bl_cnts_dnst", sep="_"),
        x_axis="blacklisted_fraction",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$maxblacklisted,
        x_label="Fraction of reads within blacklisted regions per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Density of fraction of reads within blacklisted regions per cell (", suffix, ")", sep=""),
        scale_x_log10=FALSE,
        zoom_on_intercept=TRUE,
        show_ranked=TRUE,
        pdf=args$pdf
    )
    export_geom_point_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gex_atac_umi_corr", sep="_"),
        facet_by="new.ident",
        x_axis="nCount_ATAC",
        x_label="ATAC UMIs per cell",
        y_axis="nCount_RNA",
        y_label="GEX UMIs per cell",
        x_left_intercept=args$atacminumi,
        y_low_intercept=args$gexminumi,
        alpha_intercept=1,
        color_by="mito_percentage",
        colors=c("lightslateblue", "red", "green"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        legend_title="Mitochondrial %",
        plot_title=paste("GEX vs ATAC UMIs per cell correlation (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        pdf=args$pdf
    )
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "ATAC"
    export_tss_plot(
        data=seurat_data,
        split_by="new.ident",
        group_by_value=args$mintssenrich,
        combine_guides="collect",
        rootname=paste(args$output, suffix, "tss_enrch", sep="_"),
        plot_title=paste("TSS Enrichment Score (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_fragments_hist(
        data=seurat_data,
        split_by="new.ident",
        group_by_value=args$maxnuclsignal,
        combine_guides="collect",
        rootname=paste(args$output, suffix, "frg_len_hist", sep="_"),
        plot_title=paste("Fragments Length Histogram (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    DefaultAssay(seurat_data) <- backup_assay
    export_geom_point_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gene_umi_corr", sep="_"),
        x_axis="nCount_RNA",
        y_axis="nFeature_RNA",
        facet_by="new.ident",
        x_left_intercept=args$gexminumi,
        y_low_intercept=args$mingenes,
        y_high_intercept=args$maxgenes,
        color_by="mito_percentage",
        colors=c("lightslateblue", "red", "green"),
        color_limits=c(0, 100),
        color_break=args$maxmt,
        x_label="UMIs per cell",
        y_label="Genes per cell",
        legend_title="Mitochondrial %",
        plot_title=paste("Genes vs GEX UMIs per cell correlation (", suffix, ")", sep=""),
        scale_x_log10=TRUE,
        scale_y_log10=TRUE,
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "mito_perc_dnst", sep="_"),
        x_axis="mito_percentage",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$maxmt,
        x_label="Percentage of transcripts mapped to mitochondrial genes per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Density of transcripts mapped to mitochondrial genes per cell (", suffix, ")", sep=""),
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    export_feature_scatter_plot(
        data=seurat_data,
        x_axis="nFeature_RNA",
        y_axis="mito_percentage",
        x_label="Genes",
        y_label="Mitochondrial %",
        split_by="new.ident",
        color_by="miQC.keep",
        alpha=0.1,
        palette=c("keep"="lightslateblue", "discard"="red"),
        plot_title=paste("MiQC prediction of the compromised cells level (", suffix, ")", sep=""),
        legend_title="Category",
        combine_guides="collect",
        rootname=paste(args$output, suffix, "miqc_mtrcs", sep="_"),
        pdf=args$pdf
    )
    export_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "nvlt_score_dnst", sep="_"),
        x_axis="log10_gene_per_log10_umi",
        color_by="new.ident",
        facet_by="new.ident",
        x_left_intercept=args$minnovelty,
        x_label="log10 Genes / log10 UMIs per cell",
        y_label="Density",
        legend_title="Identity",
        plot_title=paste("Novelty score density per cell (", suffix, ")", sep=""),
        zoom_on_intercept=TRUE,
        pdf=args$pdf
    )
    # export_data(
    #     get_density_data(      # if raises any exception, will be caught by export_data
    #         data=seurat_data,
    #         features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal", "frip", "blacklisted_fraction"),
    #         labels=c("GEX UMIs", "Genes", "Mitochondrial %", "Novelty score", "ATAC UMIs", "Peaks", "TSS enrichment score", "Nucleosome signal", "FRiP", "Fr. of reads in bl-ted reg."),
    #         assay="RNA",
    #         slot="data"
    #     ),
    #     location=paste(args$output, suffix, "qc_mtrcs.tsv", sep="_")
    # )
    export_vln_plot(
        data=seurat_data,
        features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment",       "nucleosome_signal", "frip", "blacklisted_fraction"),
        labels=c(  "GEX UMIs",   "Genes",        "Mitochondrial %", "Novelty score",            "ATAC UMIs",   "Peaks",         "TSS enrichment score", "Nucleosome signal", "FRiP", "Fr. of reads in bl-ted reg."),
        rootname=paste(args$output, suffix, "qc_mtrcs", sep="_"),
        plot_title=paste("QC metrics densities per cell (", suffix, ")", sep=""),
        legend_title="Identity",
        hide_x_text=TRUE,
        pt_size=0,
        palette="Paired",
        combine_guides="collect",
        pdf=args$pdf
    )
}


export_all_clustering_plots <- function(seurat_data, suffix, args) {
    cluster_prefix <- "wsnn"
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        export_dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title=paste("Clustered UMAP projected PCA of filtered GEX datasets. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "gex_umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="rnaumap",
            plot_title=paste("Split by condition clustered UMAP projected PCA of filtered GEX datasets. Resolution", current_resolution),
            legend_title="Cluster",
            split_by="condition",
            group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),
            perc_split_by="condition",
            perc_group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),            
            label=TRUE,
            rootname=paste(args$output, suffix, "gex_umap_spl_by_cond_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="atacumap",
            plot_title=paste("Clustered UMAP projected LSI of filtered ATAC datasets. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "atac_umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="atacumap",
            plot_title=paste("Split by condition clustered UMAP projected LSI of filtered ATAC datasets. Resolution", current_resolution),
            legend_title="Cluster",
            split_by="condition",
            group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),
            perc_split_by="condition",
            perc_group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),            
            label=TRUE,
            rootname=paste(args$output, suffix, "atac_umap_spl_by_cond_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="wnnumap",
            plot_title=paste("Clustered UMAP projected WNN. Resolution", current_resolution),
            legend_title="Cluster",
            group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "wnn_umap_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="wnnumap",
            plot_title=paste("Split by condition clustered UMAP projected WNN. Resolution", current_resolution),
            legend_title="Cluster",
            split_by="condition",
            group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),
            perc_split_by="condition",
            perc_group_by=paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep="."),            
            label=TRUE,
            rootname=paste(args$output, suffix, "wnn_umap_spl_by_cond_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        Idents(seurat_data) <- paste(paste(cluster_prefix, "res", sep="_"), current_resolution, sep=".")
        export_feature_plot(
            data=seurat_data,
            features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal", "frip", "blacklisted_fraction"),
            labels=c("GEX UMIs", "Genes", "Mitochondrial %", "Novelty score", "ATAC UMIs", "Peaks", "TSS enrichment score", "Nucleosome signal", "FRiP", "Fr. of reads in bl-ted reg."),
            reduction="wnnumap",
            plot_title=paste("QC metrics for clustered UMAP projected WNN. Resolution", current_resolution),
            label=TRUE,
            alpha=0.4,
            max_cutoff="q99",  # to prevent outlier cells to distort coloring
            rootname=paste(args$output, suffix, "wnn_qc_mtrcs_res", current_resolution, sep="_"),
            combine_guides="keep",
            pdf=args$pdf
        )
        Idents(seurat_data) <- "new.ident"
    }
}


export_all_dimensionality_plots <- function(seurat_data, suffix, args) {
    export_elbow_plot(
        data=seurat_data,
        ndims=50,
        rootname=paste(args$output, suffix, "gex_elbow", sep="_"),
        plot_title="Elbow plot from GEX PCA of filtered integrated/scaled datasets",
        pdf=args$pdf
    )
    export_dim_plot(
        data=seurat_data,
        reduction="pca",
        plot_title="GEX PCA of filtered integrated/scaled datasets",
        legend_title="Identity",
        rootname=paste(args$output, suffix, "gex_pca", sep="_"),
        palette="Paired",
        pdf=args$pdf
    )
    export_corr_plot(
        data=seurat_data,
        reduction="pca",
        qc_columns=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi"),
        qc_labels=c("GEX UMIs", "Genes", "Mitochondrial %", "Novelty score"),
        plot_title="Correlation plots between main QC metrics and PCA reduction on GEX assay",
        rootname=paste(args$output, suffix, "gex_qc_dim_corr", sep="_"),
        combine_guides="collect",
        pdf=args$pdf
    )
    export_corr_plot(
        data=seurat_data,
        reduction="atac_lsi",
        qc_columns=c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal", "frip", "blacklisted_fraction"),
        qc_labels=c("ATAC UMIs", "Peaks", "TSS enrichment score", "Nucleosome signal", "FRiP", "Fr. of reads in bl-ted reg."),
        plot_title="Correlation plots between main QC metrics and LSI reduction on ATAC assay",
        rootname=paste(args$output, suffix, "atac_qc_dim_corr", sep="_"),
        combine_guides="collect",
        pdf=args$pdf
    )
}


export_dot_plot <- function(data, features, rootname, plot_title, x_label, y_label, cluster_idents=FALSE, min_pct=0.01, col_min=-2.5, col_max=2.5, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plot <- DotPlot(
                        data,
                        features=features,
                        cluster.idents=cluster_idents,
                        dot.min=min_pct,
                        col.min=col_min,
                        col.max=col_max,
                        scale=TRUE,
                        scale.by="size"  # for optimal perception
                    ) +
                    xlab(x_label) +
                    ylab(y_label) +
                    theme_gray() +
                    ggtitle(plot_title) +
                    RotatedAxis()

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(plot))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(plot))
                dev.off()
            }

            print(paste("Export dot plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export dot plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_feature_plot <- function(data, features, labels, rootname, reduction, plot_title, split_by=NULL, label=FALSE, order=FALSE, min_cutoff=NA, max_cutoff=NA, pt_size=NULL, combine_guides=NULL, alpha=NULL, pdf=FALSE, width=1200, height=800, resolution=100){
    tryCatch(
        expr = {
            plots <- FeaturePlot(
                        data,
                        features=features,
                        pt.size=pt_size,
                        order=order,
                        min.cutoff=min_cutoff,
                        max.cutoff=max_cutoff,
                        reduction=reduction,
                        split.by=split_by,
                        label=label,
                        combine=FALSE       # to return a list of gglots
                    )
            plots <- lapply(seq_along(plots), function(i){
                plots[[i]] <- plots[[i]] + ggtitle(labels[i]) + theme_gray()
                if (!is.null(alpha)) { plots[[i]]$layers[[1]]$aes_params$alpha <- alpha }
                return (plots[[i]])
            })
            combined_plots <- wrap_plots(plots, guides=combine_guides) + plot_annotation(title=plot_title)

            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            suppressMessages(print(combined_plots))
            dev.off()

            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                suppressMessages(print(combined_plots))
                dev.off()
            }

            print(paste("Export feature plot to ", rootname, ".(png/pdf)", sep=""))
        },
        error = function(e){
            tryCatch(expr={dev.off()}, error=function(e){print(paste("Called  dev.off() with error -", e))})
            print(paste("Failed to export feature plot to ", rootname, ".(png/pdf)", sep=""))
        }
    )
}


export_all_expression_plots <- function(seurat_data, suffix, args, assay="RNA") {
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- assay
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        Idents(seurat_data) <- paste("wsnn_res", current_resolution, sep=".")
        export_dot_plot(
            data=seurat_data,
            features=args$gexfeatures,
            plot_title=paste("Scaled average log normalized gene expression per cluster. Resolution", current_resolution),
            x_label="Genes",
            y_label="Clusters",
            cluster_idents=FALSE,
            rootname=paste(args$output, suffix, "avg_per_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_feature_plot(
            data=seurat_data,
            features=args$gexfeatures,
            labels=args$gexfeatures,
            reduction="wnnumap",
            plot_title=paste("Log normalized gene expression per cell of clustered datasets. Resolution", current_resolution),
            label=TRUE,
            order=TRUE,
            max_cutoff="q99",  # to prevent cells with overexpressed gene from distoring the color bar
            rootname=paste(args$output, suffix, "per_clst_cell_res", current_resolution, sep="_"),
            combine_guides="keep",
            pdf=args$pdf
        )
        export_vln_plot(
            data=seurat_data,
            features=args$gexfeatures,
            labels=args$gexfeatures,
            plot_title=paste("Log normalized gene expression densities per cluster. Resolution", current_resolution),
            legend_title="Cluster",
            log=TRUE,
            pt_size=0,
            combine_guides="collect",
            rootname=paste(args$output, suffix, "dnst_per_clst_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        Idents(seurat_data) <- "orig.ident"
    }
    DefaultAssay(seurat_data) <- backup_assay
}


export_rds <- function(data, location){
    tryCatch(
        expr = {
            saveRDS(data, location)
            print(paste("Export data as RDS to", location, sep=" "))
        },
        error = function(e){
            print(paste("Failed to export data as RDS to", location, sep=" "))
        }
    )
}


load_cell_identity_data <- function(location) {
    cell_identity_data <- read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    # prepend with LETTERS, otherwise the order on the plot will be arbitrary sorted
    cell_identity_data <- cell_identity_data %>% mutate("library_id"=paste(LETTERS[1:nrow(cell_identity_data)], .$library_id))
    return (cell_identity_data)
}


load_condition_data <- function(location, cell_identity_data) {
    default_condition_data <- data.frame(
        library_id=cell_identity_data$library_id,
        condition=cell_identity_data$library_id,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    if (!is.null(location)){
        condition_data <- read.table(
            location,
            sep=get_file_type(location),
            header=TRUE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )
        # prepend with LETTERS to correspond to the library_id from the cell_identity_data
        condition_data <- condition_data %>% mutate("library_id"=paste(LETTERS[1:nrow(condition_data)], .$library_id))
        if ( (nrow(condition_data) == nrow(cell_identity_data)) && all(is.element(cell_identity_data$library_id, condition_data$library_id)) ){
            print(paste("Condition data is successfully loaded from ", location))
            return (condition_data)
        } else {
            print(paste("Applying defaults - condition data loaded from", location, "is malformed"))
            return (default_condition_data)
        }
    }
    print("Condition data is not provided. Applying defaults")
    return (default_condition_data)
}


load_barcodes_data <- function(location, seurat_data) {
    default_barcodes_data <- Cells(seurat_data)                               # to include all available cells
    if (!is.null(location)){
        barcodes_data <- read.table(
            location,
            sep=get_file_type(location),
            header=FALSE,
            check.names=FALSE,
            stringsAsFactors=FALSE
        )[,1]                                                                 # to get it as vector
        print(
            paste(
                "Barcodes data is successfully loaded and from", location,
                "and will be used to prefilter feature-barcode matrices by",
                "cells of interest")
        )
        return (barcodes_data)
    }
    print("Barcodes data is not provided. Using all cells")
    return (default_barcodes_data)
}


load_seurat_data <- function(args, cell_identity_data, condition_data) {
    raw_data <- Read10X(data.dir=args$mex) 
    seurat_data <- CreateSeuratObject(
        counts=raw_data$`Gene Expression`,
        min.cells=args$gexmincells,
        names.delim="-",
        names.field=2
    )

    annotation <- rtracklayer::import(args$annotations, format="GFF")

    seurat_data[["ATAC"]] <- CreateChromatinAssay(
        counts=raw_data$Peaks,
        sep=c(":", "-"),
        fragments=args$fragments,
        min.cells=args$atacmincells,
        annotation=annotation
    )

    if (args$callpeaks){
        print("Replacing Cell Ranger's peaks with MACS2 peaks.")
        backup_assay <- DefaultAssay(seurat_data)
        DefaultAssay(seurat_data) <- "ATAC"
        # We use group.by="orig.ident" to force CallPeaks using only those fragments
        # that belong to the cells already identified by Cell Ranger as valid.
        # Otherwise, CallPeaks function will use all fragments even from invalid cells.
        # Grouping by orig.ident will results in only 1 group, which is fine as long
        # as our input mex matrix includes only one dataset (not aggregated). In case
        # --mex pointed to aggregated datasets, for each of them MACS2 will be called
        # independently and results will be combined. This may potentially cause problem
        # when integrating multiple datasets as all the peaks should be identical between
        # them.
        macs2_peaks <- CallPeaks(seurat_data, group.by="orig.ident")
        macs2_counts <- FeatureMatrix(
            fragments=Fragments(seurat_data),
            sep=c(":", "-"),
            features=macs2_peaks,
            cells=colnames(seurat_data),
            verbose=args$verbose,
        )
        atac_assay <- CreateChromatinAssay(
            counts=macs2_counts,
            sep=c(":", "-"),
            fragments=Fragments(seurat_data),
            min.cells=args$atacmincells,
            min.features=-1,              # as they check ncount.cell > min.features and by default it's 0, we will remove cells without peaks and won't be able to add new assay to our seurat_data
            annotation=annotation
        )
        seurat_data[["ATAC"]] <- atac_assay
        DefaultAssay(seurat_data) <- backup_assay
    }

    print("Assigning new dataset identities")
    Idents(seurat_data) <- "orig.ident"                                  # safety measure to make sure we get correct Idents
    idents <- as.numeric(as.character(Idents(seurat_data)))              # need to properly convert factor to numeric vector
    new_ident <- cell_identity_data$library_id[idents]
    if (sum(is.na(new_ident)) > 0){
        print("Identity file includes less than expected number of rows. Exiting.")
        quit(save="no", status=1, runLast=FALSE)
    }
    seurat_data[["new.ident"]] <- new_ident
    seurat_data[["condition"]] <- condition_data$condition[match(seurat_data$new.ident, condition_data$library_id)]
    Idents(seurat_data) <- "new.ident"
    if ( nrow(cell_identity_data) > length(unique(as.vector(as.character(Idents(seurat_data))))) ){
        print("Identity file includes more than expected number of rows. Exiting.")
        quit(save="no", status=1, runLast=FALSE)
    }
    return (seurat_data)
}


load_blacklisted_data <- function(location) {
    default_blacklisted_data <- NULL
    if (!is.null(location)){
        blacklisted_data <- rtracklayer::import(location, format="BED")
        print(paste("Blacklisted regions data is successfully loaded from ", location))
        return (blacklisted_data)
    }
    print("Blacklisted regions data is not provided. Cells won't be filtered by --maxblacklisted")
    return (default_blacklisted_data)
}


get_args <- function(){
    parser <- ArgumentParser(description="Runs Seurat Weighted Nearest Neighbor Analysis")
    parser$add_argument(
        "--mex",
        help=paste(
            "Path to the folder with feature-barcode matrix from Cell Ranger ARC Count/Aggregate",
            "experiment in MEX format. The rows consist of all the gene and peak features",
            "concatenated together and the columns are restricted to those barcodes that are",
            "identified as cells."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--identity",
        help=paste(
            "Path to the metadata TSV/CSV file to set the datasets identities.",
            "If --mex points to the Cell Ranger ARC Aggregate outputs, the aggr.csv",
            "file can be used. If Cell Ranger ARC Count outputs have been used in",
            "--mex, the file should include at least one column - 'library_id' and",
            "one row with the alias for Cell Ranger ARC Count experiment."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--fragments",
        help=paste(
            "Count and barcode information for every ATAC fragment observed in",
            "the experiment in TSV format. Tbi-index file is required."
        ),
        type="character", required="True"
    )
    parser$add_argument(
        "--annotations",
        help="Path to the genome annotation file in GTF format",
        type="character", required="True"
    )
    parser$add_argument(
        "--condition",
        help=paste(
            "Path to the TSV/CSV file to define datasets grouping. First column -",
            "'library_id' with the values provided in the same order as in the",
            "correspondent column of the --identity file, second column 'condition'.",
            "Default: each dataset is assigned to a separate group."
        ),
        type="character"
    )
    parser$add_argument(
        "--metadata",
        help=paste(
            "Path to the TSV/CSV file to optionally extend cells metadata with categorical",
            "values using cells barcodes. First column - 'barcode' should include cells",
            "barcodes that correspond to the data provided in --mex. Values from",
            "all other columns will be added as extra metadata columns prefixed",
            "with 'custom_'. Values for missing barcodes will be set to 'Unknown'.",
            "Default: no extra cells metadata is added"
        ),
        type="character"
    )
    parser$add_argument(
        "--blacklisted",
        help="Path to the blacklisted regions file in BED format",
        type="character"
    )
    parser$add_argument(
        "--barcodes",
        help=paste(
            "Path to the headerless TSV/CSV file with the list of barcodes to select",
            "cells of interest (one barcode per line). Prefilters input feature-barcode",
            "matrix to include only selected cells.",
            "Default: use all cells."
        ),
        type="character"
    )
    parser$add_argument(
        "--gexmincells",
        help=paste(
            "Include only GEX features detected in at least this many cells.",
            "Default: 5 (applied to all datasets)"
        ),
        type="integer", default=5
    )
    parser$add_argument(
        "--mingenes",
        help=paste(
            "Include cells where at least this many GEX features are detected.",
            "If multiple values provided, each of them will be applied to the",
            "correspondent dataset from the --mex input based on the --identity",
            "file.",
            "Default: 250 (applied to all datasets)"
        ),
        type="integer", default=250, nargs="*"
    )
    parser$add_argument(
        "--maxgenes",
        help=paste(
            "Include cells with the number of GEX features not bigger than this value.",
            "If multiple values provided, each of them will be applied to the correspondent",
            "dataset from the --mex input based on the --identity file.",
            "Default: 5000 (applied to all datasets)"
        ),
        type="integer", default=5000, nargs="*"
    )
    parser$add_argument(
        "--gexminumi",
        help=paste(
            "Include cells where at least this many GEX UMIs (transcripts) are detected.",
            "If multiple values provided, each of them will be applied to the correspondent",
            "dataset from the --mex input based on the --identity file.",
            "Default: 500 (applied to all datasets)"
        ),
        type="integer", default=500, nargs="*"
    )
    parser$add_argument(
        "--mitopattern",
        help=paste(
            "Regex pattern to identify mitochondrial GEX features.",
            "Default: '^Mt-'"
        ),
        type="character", default="^Mt-"
    )
    parser$add_argument(
        "--maxmt",
        help=paste(
            "Include cells with the percentage of GEX transcripts mapped to mitochondrial",
            "genes not bigger than this value.",
            "Default: 5 (applied to all datasets)"
        ),
        type="double", default=5
    )
    parser$add_argument(
        "--regressmt",
        help=paste(
            "Regress mitochondrial genes expression as a confounding source of variation.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--minnovelty",
        help=paste(
            "Include cells with the novelty score not lower than this value, calculated for",
            "GEX as log10(genes)/log10(UMIs). If multiple values provided, each of them will",
            "be applied to the correspondent dataset from the --mex input based on the",
            "--identity file.",
            "Default: 0.8 (applied to all datasets)"
        ),
        type="double", default=0.8, nargs="*"
    )
    parser$add_argument(
        "--atacmincells",
        help=paste(
            "Include only ATAC features detected in at least this many cells.",
            "Default: 5 (applied to all datasets)"
        ),
        type="integer", default=5
    )
    parser$add_argument(
        "--atacminumi",
        help=paste(
            "Include cells where at least this many ATAC UMIs (transcripts) are detected.",
            "If multiple values provided, each of them will be applied to the correspondent",
            "dataset from the --mex input based on the --identity file.",
            "Default: 1000 (applied to all datasets)"
        ),
        type="integer", default=1000, nargs="*"
    )
    parser$add_argument(
        "--maxnuclsignal",
        help=paste(
            "Include cells with the nucleosome signal not bigger than this value.",
            "Nucleosome signal quantifies the approximate ratio of mononucleosomal",
            "to nucleosome-free fragments. If multiple values provided, each of",
            "them will be applied to the correspondent dataset from the --mex input",
            "based on the --identity file",
            "Default: 4 (applied to all datasets)"
        ),
        type="double", default=4, nargs="*"
    )
    parser$add_argument(
        "--mintssenrich",
        help=paste(
            "Include cells with the TSS enrichment score not lower than this value.",
            "Score is calculated based on the ratio of fragments centered at the TSS",
            "to fragments in TSS-flanking regions. If multiple values provided, each",
            "of them will be applied to the correspondent dataset from the --mex input",
            "based on the --identity file.",
            "Default: 2 (applied to all datasets)"
        ),
        type="double", default=2, nargs="*"
    )
    parser$add_argument(
        "--minfrip",
        help=paste(
            "Include cells with the FRiP not lower than this value. If multiple values",
            "provided, each of them will be applied to the correspondent dataset from",
            "the --mex input based on the --identity file.",
            "Default: 0.15 (applied to all datasets)"
        ),
        type="double", default=0.15, nargs="*"
    )
    parser$add_argument(
        "--maxblacklisted",
        help=paste(
            "Include cells with the ratio of reads in genomic blacklist regions",
            "not bigger than this value. If multiple values provided, each of them",
            "will be applied to the correspondent dataset from the --mex input based",
            "on the --identity file.",
            "Default: 0.05 (applied to all datasets)"
        ),
        type="double", default=0.05, nargs="*"
    )
    parser$add_argument(
        "--callpeaks",
        help=paste(
            "Call peaks with MACS2 instead of those that are provided by Cell Ranger ARC Count.",
            "If --mex points to the Cell Ranger ARC Aggregate experiment, peaks will be called for",
            "each dataset independently and then combined",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--gexfeatures",
        help=paste(
            "GEX features of interest to evaluate expression.",
            "Default: None"
        ),
        type="character", nargs="*"
    )
    parser$add_argument(
        "--highvargex",
        help=paste(
            "Number of highly variable GEX features to detect. Used for GEX datasets",
            "integration, scaling, and dimensional reduction.",
            "Default: 3000"
        ),
        type="integer", default=3000
    )
    parser$add_argument(
        "--gexndim",
        help=paste(
            "Dimensionality to use in GEX UMAP projection and clustering (from 1 to 50).",
            "If single number N is provided, use from 1 to N PCs. If multiple numbers are",
            "provided, subset to only selected PCs.",
            "Default: from 1 to 10"
        ),
        type="integer", default=10, nargs="*"
    )
    parser$add_argument(
        "--gexlogfc",
        help=paste(
            "For putative gene markers identification include only those GEX features that",
            "on average have log fold change difference in expression between every tested",
            "pair of clusters not lower than this value.",
            "Default: 0.25"
        ),
        type="double", default=0.25
    )
    parser$add_argument(
        "--gexminpct",
        help=paste(
            "For putative gene markers identification include only those GEX features that",
            "are detected in not lower than this fraction of cells in either of the two",
            "tested clusters.",
            "Default: 0.1"
        ),
        type="double", default=0.1
    )
    parser$add_argument(
        "--gexonlypos",
        help=paste(
            "For putative gene markers identification return only positive markers.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--gextestuse",
        help=paste(
            "Statistical test to use for putative gene markers identification.",
            "Default: wilcox"
        ),
        type="character", default="wilcox",
        choices=c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
    )
    parser$add_argument(
        "--nosct",
        help=paste(
            "Do not use SCTransform when running RNA datasets integration.",
            "Use LogNormalize instead. Ignored when --mex points to the",
            "Cell Ranger ARC Count outputs (single, not aggregated dataset",
            "that doesn't require any integration) or --skipgexntrg parameter",
            "was applied.",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--atacndim",
        help=paste(
            "Dimensionality to use in ATAC UMAP projection and clustering (from 2 to 50).",
            "If single number N is provided, use from 2 to N LSI components. If multiple",
            "numbers are provided, subset to only selected LSI components.",
            "Default: from 2 to 10"
        ),
        type="integer", default=10, nargs="*"
    )
    parser$add_argument(
        "--highvaratac",
        help=paste(
            "Minimum percentile to set the top most common ATAC features as highly variable.",
            "For example, setting to 5 will use the the top 95 percent most common among all cells",
            "ATAC features as highly variable. Used for ATAC datasets integration, scaling,",
            "and dimensional reduction.",
            "Default: 75 (use only the top 25 percent of all common peaks)"
        ),
        type="integer", default=75
    )
    parser$add_argument(
        "--resolution",
        help=paste(
            "Clustering resolution. Can be set as an array.",
            "Default: 0.3, 0.5, 1.0"
        ),
        type="double", default=c(0.3, 0.5, 1.0), nargs="*"
    )
    parser$add_argument(
        "--skipgexntrg",
        help=paste(
            "Do not integrate RNA datasets, use merged data instead.",
            "Applied by default if --mex points to the Cell Ranger ARC Count",
            "outputs (single, not aggregated dataset that doesn't require any",
            "integration).",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--skipatacntrg",
        help=paste(
            "Do not integrate ATAC datasets, use merged data instead.",
            "Applied by default if --mex pointed to the Cell Ranger ARC Count",
            "outputs (single, not aggregated dataset that doesn't require any",
            "integration).",
            "Default: false"
        ),
        action="store_true"
    )
    # Export results
    parser$add_argument(
        "--pdf",
        help="Export plots in PDF. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--rds",
        help="Save Seurat data to RDS file. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--verbose",
        help="Print debug information. Default: false",
        action="store_true"
    )
    parser$add_argument(
        "--skipmiqc",
        help=paste(
            "Skip threshold prediction for the percentage of transcripts",
            "mapped to mitochondrial genes (do not run MiQC).",
            "Default: false"
        ),
        action="store_true"
    )
    parser$add_argument(
        "--output",
        help="Output prefix. Default: ./seurat",
        type="character", default="./seurat"
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
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}

args <- get_args()

cat("Step 0: Runtime configuration\n")
print(args)

print(
    paste(
        "Setting parallelization parameters to", args$cpus,
        "cores, and", args$memory, "GB of memory")
    )
setup_parallelization(args)

cat("\n\nStep 1: Loading raw datasets\n")

print(paste("Loading datasets identity data from", args$identity))
cell_identity_data <- load_cell_identity_data(args$identity)

print("Trying to load condition data")
condition_data <- load_condition_data(args$condition, cell_identity_data)

print(paste("Loading gene/peak-barcode matrices from", args$mex))
print(paste("Loading fragments from", args$fragments))
print(paste("Loading annotations from from", args$annotations))
seurat_data <- load_seurat_data(args, cell_identity_data, condition_data)
debug_seurat_data(seurat_data, args)

print("Validating/adjusting input parameters")
idents_count <- length(unique(as.vector(as.character(Idents(seurat_data)))))
for (key in names(args)){
    if (key %in% c("mingenes", "maxgenes", "gexminumi", "minnovelty", "atacminumi", "maxnuclsignal", "mintssenrich", "minfrip", "maxblacklisted")){
        if (length(args[[key]]) != 1 && length(args[[key]]) != idents_count){
            print(paste("Filtering parameter", key, "has an ambiguous size. Exiting"))
            quit(save="no", status=1, runLast=FALSE)
        }
        if (length(args[[key]]) == 1){
            print(paste("Extending filtering parameter", key, "to have a proper size"))
            args[[key]] <- rep(args[[key]][1], idents_count)
        }
    }
}
args$highvaratac <- paste0("q", args$highvaratac)          # need to have it in a form of "q75", for example
if (length(args$gexndim) == 1) {                           # only one value was provided, so we need to inflate it to 1:N
    args$gexndim <- c(1:args$gexndim[1])
}
if (length(args$atacndim) == 1) {                          # only one value was provided, so we need to inflate it to 2:N
    args$atacndim <- c(2:args$atacndim[1])                 # skip first LSI component by default
}
print("Adjusted parameters")
print(args)

print("Trying to load blacklisted regions")
blacklisted_data <- load_blacklisted_data(args$blacklisted)

print("Trying to load barcodes to prefilter feature-barcode matrices by cells of interest")
barcodes_data <- load_barcodes_data(args$barcodes, seurat_data)
seurat_data <- apply_cell_filters(seurat_data, barcodes_data)

print("Trying to extend cells metadata with categorical values using cells barcodes")
seurat_data <- extend_metadata_by_barcode(seurat_data, args$metadata)

print("Adding QC metrics to raw seurat data")
seurat_data <- add_qc_metrics(seurat_data, blacklisted_data, args)
debug_seurat_data(seurat_data, args)
export_all_qc_plots(seurat_data, "raw", args)

cat("\n\nStep 2: Filtering raw datasets\n")
seurat_data <- apply_qc_filters(seurat_data, cell_identity_data, args)
debug_seurat_data(seurat_data, args)
export_all_qc_plots(seurat_data, "fltr", args)

cat("\n\nStep 3: Running GEX analysis\n")
print("Running datasets integration/scaling on RNA assay")
seurat_data <- integrate_gex_data(seurat_data, args)                                            # sets "gex_integrated" as a default assay for integrated data, and leave "RNA" as a default assays for scaled data
debug_seurat_data(seurat_data, args)
print(
    paste(
        "Performing PCA reduction on", DefaultAssay(seurat_data),
        "assay using 50 principal components"
    )
)
seurat_data <- RunPCA(seurat_data, npcs=50, verbose=args$verbose)
seurat_data <- RunUMAP(
    seurat_data,
    reduction="pca",
    dims=args$gexndim,
    reduction.name="rnaumap",
    reduction.key="RNAUMAP_",
    verbose=args$verbose
)
debug_seurat_data(seurat_data, args)

cat("\n\nStep 4: Running ATAC analysis\n")
print("Running datasets integration/scaling on ATAC assay")
seurat_data <- integrate_atac_data(seurat_data, args)                # shouldn't change default assay, but will add "atac_lsi" reduction
debug_seurat_data(seurat_data, args)
seurat_data <- RunUMAP(
    seurat_data,
    reduction="atac_lsi",
    dims=args$atacndim,
    reduction.name="atacumap",
    reduction.key="ATACUMAP_",
    verbose=args$verbose
)
debug_seurat_data(seurat_data, args)
export_all_dimensionality_plots(seurat_data, "ntgr", args)

cat("\n\nStep 5: Running WNN analysis\n")
seurat_data <- FindMultiModalNeighbors(
    seurat_data,
    reduction.list=list("pca", "atac_lsi"),          # we used "atac_lsi" reduction name for both after integration and just scaling/normalization
    dims.list=list(args$gexndim, args$atacndim),
    snn.graph.name="wsnn",
    weighted.nn.name="weighted.nn",
    verbose=args$verbose
)
seurat_data <- RunUMAP(
    seurat_data,
    nn.name="weighted.nn",
    reduction.name="wnnumap",
    reduction.key="WNNUMAP_",
    verbose=args$verbose
)
seurat_data <- FindClusters(
    seurat_data,
    graph.name="wsnn",
    algorithm=3,
    resolution=args$resolution,
    verbose=args$verbose
)
debug_seurat_data(seurat_data, args)
export_all_clustering_plots(seurat_data, "clst", args)

if ("gex_integrated" %in% names(seurat_data@assays)) {                                                         # if we run integration, our counts and data slots in RNA assay remain the same, so we need to normalize counts first
    print("Normalizing counts in RNA assay")                                                                   # in case is was already normalized, it's safe to run it again, as it uses counts and overwrites data slots
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "RNA"
    seurat_data <- NormalizeData(seurat_data, verbose=args$verbose)
    DefaultAssay(seurat_data) <- backup_assay
}

if (args$rds){
    DefaultAssay(seurat_data) <- "RNA"
    export_rds(seurat_data, paste(args$output, "_clst_data.rds", sep=""))
}

if (!is.null(args$gexfeatures)){
    print("Check genes of interest to include only those that are present in the datasets")
    backup_assay <- DefaultAssay(seurat_data)
    DefaultAssay(seurat_data) <- "RNA"
    args$gexfeatures <- unique(args$gexfeatures)
    args$gexfeatures <- args$gexfeatures[args$gexfeatures %in% as.vector(as.character(rownames(seurat_data)))]
    print(args$gexfeatures)
    DefaultAssay(seurat_data) <- backup_assay
}

if (!is.null(args$gexfeatures)){ export_all_expression_plots(seurat_data, "expr", args, assay="RNA") }

print("Identifying putative gene markers for all clusters and all resolutions")
all_putative_markers <- get_all_putative_markers(seurat_data, args)
export_data(all_putative_markers, paste(args$output, "_clst_pttv_gene_markers.tsv", sep=""))

print("Exporting UCSC Cellbrowser data")
export_cellbrowser_data(
    seurat_data=seurat_data,
    assay="RNA",
    matrix_slot="counts",
    resolution=args$resolution,
    features=args$gexfeatures,
    rootname=paste(args$output, "_cellbrowser", sep="")
)