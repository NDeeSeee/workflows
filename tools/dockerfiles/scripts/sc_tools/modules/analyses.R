import("dplyr", attach=FALSE)
import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)
import("tibble", attach=FALSE)
import("glmGamPoi", attach=FALSE)  # safety measure. we don't use it directly, but SCTransform with method="glmGamPoi" needs it
import("magrittr", `%>%`, attach=TRUE)

export(
    "gex_analyze",
    "gex_cluster",
    "gex_preprocess",
    "gex_log_single",
    "gex_sct_single",
    "gex_log_integrated",
    "gex_sct_integrated",
    "get_vars_to_regress",
    "get_cell_cycle_scores",
    "gex_integration_k_weight",
    "gex_putative_gene_markers"
)

get_vars_to_regress <- function(seurat_data, args, exclude_columns=NULL) {
    vars_to_regress <- NULL
    arguments <- c(args$regressmt, args$regressgexumi,  args$regressgenes, args$regresscellcycle)   # any of then can be also NULL
    metadata_columns <- c("mito_percentage", "nCount_RNA", "nFeature_RNA", "S.Score&G2M.Score")
    for (i in 1:length(arguments)) {
        current_argument <- arguments[i]
        current_column <- metadata_columns[i]
        if ( is.null(current_argument) || (current_column %in% exclude_columns) ) {
            next
        }
        current_column <- base::unlist(base::strsplit(metadata_columns[i], "&"))
        if ( !all(current_column %in% base::colnames(seurat_data@meta.data)) ){                            # the column doesn't exists in metadata
            next
        }
        if (current_argument) {
            if (is.null(vars_to_regress)) {
                vars_to_regress <- current_column
            } else {
                vars_to_regress <- base::append(vars_to_regress, current_column)
            }
        }
    }
    return (vars_to_regress)
}

get_cell_cycle_scores <- function(seurat_data, assay, cell_cycle_data){   # we need this function to fail if something went wrong, so no tryCatch inside
    SeuratObject::DefaultAssay(seurat_data) <- assay                      # safety measure
    seurat_data <- Seurat::CellCycleScoring(
        seurat_data,
        s.features=base::as.vector(cell_cycle_data[base::tolower(cell_cycle_data$phase)=="s", "gene_id"]),
        g2m.features=base::as.vector(cell_cycle_data[base::tolower(cell_cycle_data$phase)=="g2/m", "gene_id"]),
        assay=assay,
        verbose=FALSE
    )
    return (seurat_data)
}

gex_log_single <- function(seurat_data, args, cell_cycle_data=NULL){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                                # safety measure
    base::print("Applying LogNormalize")
    scaled_norm_seurat_data <- Seurat::NormalizeData(seurat_data, verbose=FALSE)
    if(!is.null(cell_cycle_data)){
        base::tryCatch(
            expr = {
                base::print("Trying to assign cell cycle scores for RNA assay")
                scaled_norm_seurat_data <- get_cell_cycle_scores(                                       # if succeded adds S.Score and G2M.Score columns to medadata
                    scaled_norm_seurat_data,
                    "RNA",
                    cell_cycle_data
                )
            },
            error = function(e){
                base::print(base::paste("Failed to run cell cycle scoring for RNA assay with error - ", e))
            }
        )
    }
    scaled_norm_seurat_data <- Seurat::FindVariableFeatures(
        scaled_norm_seurat_data,
        nfeatures=args$highvargex,
        verbose=FALSE
    )
    vars_to_regress <- get_vars_to_regress(scaled_norm_seurat_data, args)                           # may or may not include S.Score and G2M.Score columns
    base::print(base::paste("Regressing out", paste(vars_to_regress, collapse=", ")))
    scaled_norm_seurat_data <- Seurat::ScaleData(
        scaled_norm_seurat_data,
        vars.to.regress=vars_to_regress,
        verbose=FALSE
    )
    SeuratObject::DefaultAssay(scaled_norm_seurat_data) <- "RNA"
    return (scaled_norm_seurat_data)
}

gex_sct_single <- function(seurat_data, args, cell_cycle_data=NULL){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                # safety measure
    method <- base::ifelse(args$gexnorm=="sctglm", "glmGamPoi", "poisson")
    vars_to_regress <- get_vars_to_regress(seurat_data, args)                       # will never include "S.Score", "G2M.Score" columns as cell cycle scores are not assigned
    base::print(
        base::paste(
            "Applying SCTransform using", method,
            "method for initial parameter estimation.",
            "Regressing out", paste(vars_to_regress, collapse=", ")
        )
    )
    scaled_norm_seurat_data <- Seurat::SCTransform(
        seurat_data,                                                                # use not splitted seurat_data
        assay="RNA",
        new.assay.name="SCT",
        variable.features.n=args$highvargex,
        method=method,
        vars.to.regress=vars_to_regress,                                            # first portion of variables to regress. will never include "S.Score" and "G2M.Score"
        conserve.memory=args$lowmem,
        verbose=FALSE                                                               # too many stdout
    )
    if(!is.null(cell_cycle_data)){
        base::tryCatch(
            expr = {
                base::print("Trying to assign cell cycle scores for SCT assay")
                scaled_norm_seurat_data <- get_cell_cycle_scores(                       # if succeded adds S.Score and G2M.Score columns to medadata
                    scaled_norm_seurat_data,
                    "SCT",
                    cell_cycle_data
                )
                vars_to_regress <- get_vars_to_regress(scaled_norm_seurat_data, args)   # will include S.Score and G2M.Score
                base::print(
                    base::paste(
                        "Re-applying SCTransform using", method,
                        "method for initial parameter estimation.",
                        "Regressing out", paste(vars_to_regress, collapse=", ")
                    )
                )
                scaled_norm_seurat_data <- Seurat::SCTransform(
                    seurat_data,
                    assay="RNA",
                    new.assay.name="SCT",
                    variable.features.n=args$highvargex,
                    method=method,
                    vars.to.regress=vars_to_regress,
                    conserve.memory=args$lowmem,
                    verbose=FALSE
                )
            },
            error = function(e){
                base::print(base::paste("Failed to run cell cycle scoring for SCT assay with error - ", e))
            }
        )
    }
    SeuratObject::DefaultAssay(scaled_norm_seurat_data) <- "SCT"
    return (scaled_norm_seurat_data)
}

gex_integration_k_weight <- function(splitted_seurat_data, default_k_weight=100){
    # For Seurat::IntegrateData k.weight 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering

    k_weight <- min(
        min(
            table(
                unlist(
                    lapply(
                        lapply(
                            lapply(
                                splitted_seurat_data,
                                SeuratObject::Idents
                            ),
                            as.character
                        ),
                        base::as.vector
                    )
                )
            )
        ),
        default_k_weight
    )
    return (k_weight)
}


gex_log_integrated <- function(splitted_seurat_data, args, cell_cycle_data=NULL){
    failed_cell_cycle_scoring <- FALSE
    for (i in 1:length(splitted_seurat_data)){
        SeuratObject::DefaultAssay(splitted_seurat_data[[i]]) <- "RNA"                            # safety measure
        base::print(
            base::paste(
                "Applying LogNormalize for", SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                "dataset"
            )
        )
        splitted_seurat_data[[i]] <- Seurat::NormalizeData(
            splitted_seurat_data[[i]],
            verbose=FALSE
        )
        if(!is.null(cell_cycle_data)){
            if (!failed_cell_cycle_scoring){
                base::tryCatch(
                    expr = {
                        base::print(
                            base::paste(
                                "Trying to assign cell cycle scores for RNA assay of",
                                SeuratObject::Idents(splitted_seurat_data[[i]])[1], "dataset"
                            )
                        )
                        splitted_seurat_data[[i]] <- get_cell_cycle_scores(                       # if succeded adds S.Score and G2M.Score columns to medadata
                            splitted_seurat_data[[i]],
                            "RNA",
                            cell_cycle_data
                        )
                    },
                    error = function(e){
                        base::print(base::paste("Failed to run cell cycle scoring for RNA assay with error - ", e))
                        failed_cell_cycle_scoring <- TRUE
                    }
                )
            } else {
                base::print(
                    base::paste(
                        "Skip cell cycle scores assignment for RNA assay of",
                        SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                        "dataset due to failure in one of the other datasets."
                    )
                )
            }
        }
        splitted_seurat_data[[i]] <- Seurat::FindVariableFeatures(
            splitted_seurat_data[[i]],
            nfeatures=args$highvargex,
            verbose=FALSE
        )
    }
    integration_features <- Seurat::SelectIntegrationFeatures(
        splitted_seurat_data,
        nfeatures=args$highvargex,
        verbose=FALSE
    )
    integration_anchors <- Seurat::FindIntegrationAnchors(
        splitted_seurat_data,
        normalization.method="LogNormalize",
        anchor.features=integration_features,
        verbose=FALSE
    )
    integrated_seurat_data <- Seurat::IntegrateData(
        integration_anchors, 
        new.assay.name="gex_integrated",
        normalization.method="LogNormalize",
        k.weight=gex_integration_k_weight(splitted_seurat_data, 100),
        verbose=FALSE
    )
    if (failed_cell_cycle_scoring){
        base::print(
            base::paste(
                "At least one of the datasets failed in cell cycle score assignment.",
                "Removing S.Score and G2M.Score columns from metadata."
            )
        )
        integrated_seurat_data[["S.Score"]] <- NULL
        integrated_seurat_data[["G2M.Score"]] <- NULL
    }
    vars_to_regress <- get_vars_to_regress(integrated_seurat_data, args)                # may or may not include S.Score and G2M.Score columns
    base::print(base::paste("Regressing out", paste(vars_to_regress, collapse=", ")))
    integrated_seurat_data <- Seurat::ScaleData(
        integrated_seurat_data,
        vars.to.regress=vars_to_regress,
        verbose=FALSE
    )
    SeuratObject::DefaultAssay(integrated_seurat_data) <- "gex_integrated"
    base::rm(integration_features, integration_anchors)                                 # remove unused data
    base::gc(verbose=FALSE)
    return (integrated_seurat_data)
}

gex_sct_integrated <- function(splitted_seurat_data, args, cell_cycle_data=NULL){
    method <- base::ifelse(args$gexnorm=="sctglm", "glmGamPoi", "poisson")
    failed_cell_cycle_scoring <- FALSE
    for (i in 1:length(splitted_seurat_data)) {
        SeuratObject::DefaultAssay(splitted_seurat_data[[i]]) <- "RNA"                            # safety measure
        vars_to_regress <- get_vars_to_regress(splitted_seurat_data[[i]], args)                   # will never include S.Score and G2M.Score
        base::print(
            base::paste(
                "Applying SCTransform for", SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                "dataset using", method, "method for initial parameter estimation.",
                "Regressing out", paste(vars_to_regress, collapse=", ")
            )
        )
        splitted_seurat_data[[i]] <- Seurat::SCTransform(
            splitted_seurat_data[[i]],
            assay="RNA",
            new.assay.name="SCT",
            variable.features.n=args$highvargex,
            method=method,
            vars.to.regress=vars_to_regress,
            conserve.memory=args$lowmem,
            verbose=FALSE
        )
        if(!is.null(cell_cycle_data)){
            if (!failed_cell_cycle_scoring){
                base::tryCatch(
                    expr = {
                        base::print(
                            base::paste(
                                "Trying to assign cell cycle scores for SCT assay of",
                                SeuratObject::Idents(splitted_seurat_data[[i]])[1], "dataset"
                            )
                        )
                        splitted_seurat_data[[i]] <- get_cell_cycle_scores(                              # if succeded adds S.Score and G2M.Score columns to medadata
                            splitted_seurat_data[[i]],
                            "SCT",
                            cell_cycle_data
                        )
                        vars_to_regress <- get_vars_to_regress(splitted_seurat_data[[i]], args)          # will include S.Score and G2M.Score
                        base::print(
                            base::paste(
                                "Re-applying SCTransform for", SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                                "dataset using", method, "method for initial parameter estimation.",
                                "Regressing out", paste(vars_to_regress, collapse=", ")
                            )
                        )
                        splitted_seurat_data[[i]] <- Seurat::SCTransform(
                            splitted_seurat_data[[i]],
                            assay="RNA",
                            new.assay.name="SCT",
                            variable.features.n=args$highvargex,
                            method=method,
                            vars.to.regress=vars_to_regress,
                            conserve.memory=args$lowmem,
                            verbose=FALSE
                        )
                    },
                    error = function(e){
                        base::print(base::paste("Failed to run cell cycle scoring for SCT assay with error - ", e))
                        failed_cell_cycle_scoring <- TRUE
                    }
                )
            } else {
                base::print(
                    base::paste(
                        "Skip cell cycle scores assignment for SCT assay of",
                        SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                        "dataset due to failure in one of the other datasets."
                    )
                )
            }
        }
    }
    if (failed_cell_cycle_scoring){
        base::print(base::paste("At least one of the datasets failed in cell cycle score assignment."))
        for (i in 1:length(splitted_seurat_data)){
            if(!is.null(args$regresscellcycle) && args$regresscellcycle){
                vars_to_regress <- get_vars_to_regress(splitted_seurat_data[[i]], args, "S.Score&G2M.Score")    # force to exclude "S.Score&G2M.Score"
                base::print(
                    base::paste(
                        "Re-applying SCTransform for", SeuratObject::Idents(splitted_seurat_data[[i]])[1],
                        "dataset using", method, "method for initial parameter estimation.",
                        "Regressing out", paste(vars_to_regress, collapse=", "),
                        "Skip cell cycle score assignment."
                    )
                )
                splitted_seurat_data[[i]] <- Seurat::SCTransform(
                    splitted_seurat_data[[i]],
                    assay="RNA",
                    new.assay.name="SCT",
                    variable.features.n=args$highvargex,
                    method=method,
                    vars.to.regress=vars_to_regress,
                    conserve.memory=args$lowmem,
                    verbose=FALSE
                )
            } else {
                base::print(
                    base::paste(
                        "Removing S.Score and G2M.Score columns from the metadata of",
                        SeuratObject::Idents(splitted_seurat_data[[i]])[1]
                    )
                )
                splitted_seurat_data[[i]][["S.Score"]] <- NULL
                splitted_seurat_data[[i]][["G2M.Score"]] <- NULL
            }
        }
    }
    integration_features <- Seurat::SelectIntegrationFeatures(
        splitted_seurat_data,
        nfeatures=args$highvargex,
        verbose=FALSE
    )
    splitted_seurat_data <- Seurat::PrepSCTIntegration(
        splitted_seurat_data, 
        anchor.features=integration_features,
        verbose=FALSE
    )
    integration_anchors <- Seurat::FindIntegrationAnchors(
        splitted_seurat_data,
        normalization.method="SCT",
        anchor.features=integration_features,
        verbose=FALSE
    )
    integrated_seurat_data <- Seurat::IntegrateData(
        integration_anchors, 
        new.assay.name="gex_integrated",
        normalization.method="SCT",
        k.weight=gex_integration_k_weight(splitted_seurat_data, 100),
        verbose=FALSE
    )
    SeuratObject::DefaultAssay(integrated_seurat_data) <- "gex_integrated"
    base::rm(integration_features, integration_anchors)                         # remove unused data
    base::gc(verbose=FALSE)
    return (integrated_seurat_data)
}

gex_preprocess <- function(seurat_data, args, cell_cycle_data=NULL) {
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                                        # safety measure
    splitted_seurat_data <- Seurat::SplitObject(seurat_data, split.by="new.ident")                          # to check if we have aggregated datasets
    if (args$ntgr == "none" | length(splitted_seurat_data) == 1){
        base::print(
            base::paste(
                "Skip datasets integration (either forced or only one identity is present)",
                "using the original not splitted seurat data."
            )
        )
        if (args$gexnorm == "log"){
            processed_seurat_data <- gex_log_single(seurat_data, args, cell_cycle_data)                     # sets default assay to RNA
        } else {
            processed_seurat_data <- gex_sct_single(seurat_data, args, cell_cycle_data)                     # sets default assay to SCT
        }
    } else {
        base::print("Run datasets integration using splitted seurat data.")
        if (args$gexnorm == "log"){
            processed_seurat_data <- gex_log_integrated(splitted_seurat_data, args, cell_cycle_data)        # sets default assay to gex_integrated
        } else {
            processed_seurat_data <- gex_sct_integrated(splitted_seurat_data, args, cell_cycle_data)        # sets default assay to gex_integrated
        }
    }
    base::rm(splitted_seurat_data)
    base::gc(verbose=FALSE)
    return (processed_seurat_data)
}

gex_analyze <- function(seurat_data, args, cell_cycle_data=NULL){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    backup_reductions <- c()                                                    # GEX integration main remove atac related reductions so we need to back them up
    for (reduction_name in c("atac_lsi", "atacumap")){
        if (reduction_name %in% names(seurat_data@reductions)){
            base::print(base::paste("Backing up reduction", reduction_name))
            backup_reductions[[reduction_name]] <- seurat_data[[reduction_name]]
        }
    }
    seurat_data <- gex_preprocess(seurat_data, args, cell_cycle_data)           # sets "gex_integrated" as a default assay for integrated data, and either "RNA" or "SCT" for not integrated data
    if (length(backup_reductions) > 0){                                          # restoring backed up reductions
        for (reduction_name in names(backup_reductions)){
            base::print(base::paste("Restoring reduction", reduction_name, "from backup"))
            seurat_data[[reduction_name]] <- backup_reductions[[reduction_name]]
        }
    }
    base::print(
        base::paste(
            "Performing PCA reduction on", SeuratObject::DefaultAssay(seurat_data),
            "assay using 50 principal components"
        )
    )
    seurat_data <- Seurat::RunPCA(seurat_data, npcs=50, verbose=FALSE)          # add "pca" reduction that should be used in UMAP
    seurat_data <- Seurat::RunUMAP(
        seurat_data,
        reduction="pca",
        dims=args$gexndim,
        reduction.name="rnaumap",
        reduction.key="RNAUMAP_",
        spread=base::ifelse(is.null(args$uspread), 1, args$uspread),
        min.dist=base::ifelse(is.null(args$umindist), 0.3, args$umindist),
        n.neighbors=base::ifelse(is.null(args$uneighbors), 30, args$uneighbors),
        metric=base::ifelse(is.null(args$umetric), "cosine", args$umetric),
        umap.method=base::ifelse(is.null(args$umethod), "uwot", args$umethod),
        verbose=FALSE
    )
    return (seurat_data)
}

gex_cluster <- function(seurat_data, graph_name, args){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    seurat_data <- Seurat::FindNeighbors(
        seurat_data,
        annoy.metric=base::ifelse(is.null(args$ametric), "euclidean", args$ametric),
        reduction="pca",
        dims=args$gexndim,
        graph.name=base::paste(graph_name, c("_nn", ""), sep=""),
        verbose=FALSE
    )
    seurat_data <- Seurat::FindClusters(
        seurat_data,
        resolution=args$resolution,
        graph.name=graph_name,
        verbose=FALSE
    )
    return (seurat_data)
}

gex_putative_gene_markers <- function(seurat_data, resolution_prefix, args, min_diff_pct=-Inf){
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                            # safety measure
    SeuratObject::Idents(seurat_data) <- "new.ident"                            # safety measure
    all_putative_markers <- NULL
    for (i in 1:length(args$resolution)) {
        resolution <- args$resolution[i]
        base::tryCatch(
            expr = {
                SeuratObject::Idents(seurat_data) <- paste(resolution_prefix, resolution, sep=".")
                markers <- Seurat::FindAllMarkers(
                    seurat_data,
                    logfc.threshold=base::ifelse(is.null(args$gexlogfc), 0.25, args$gexlogfc),
                    min.pct=base::ifelse(is.null(args$gexminpct), 0.1, args$gexminpct),
                    only.pos=base::ifelse(is.null(args$gexonlypos), FALSE, args$gexonlypos),
                    test.use=base::ifelse(is.null(args$gextestuse), "wilcox", args$gextestuse),
                    min.diff.pct=min_diff_pct,
                    verbose=FALSE
                ) %>% dplyr::relocate(cluster, gene, .before=1)
                if (base::nrow(markers) > 0) {
                    markers <- markers %>% base::cbind(resolution=resolution, .)
                } else {
                    markers <- markers %>% tibble::add_column(resolution=base::numeric(), .before=1)  # safety measure in case markers was empty
                }
                if (!is.null(all_putative_markers)) {
                    all_putative_markers <- base::rbind(all_putative_markers, markers)
                } else {
                    all_putative_markers <- markers
                }
                base::rm(markers)                                                             # remove unused data
            },
            error = function(e){
                base::print(base::paste("Failed to identify putative gene markers for resolution", resolution, "with error -", e))
            },
            finally = {
                SeuratObject::Idents(seurat_data) <- "new.ident"
            }
        )
    }
    base::gc(verbose=FALSE)
    return (all_putative_markers)
}