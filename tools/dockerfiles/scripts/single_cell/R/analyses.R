import("Seurat", attach=FALSE)
import("Signac", attach=FALSE)

export(
    "gex_integrate",
    "gex_analyze"
)


gex_integrate <- function(seurat_data, args) {
    SeuratObject::DefaultAssay(seurat_data) <- "RNA"                                         # safety measure, in case RNA was not default assay
    splitted_seurat_data <- Seurat::SplitObject(seurat_data, split.by="new.ident")
    vars_to_regress <- NULL                                                                  # used in all three places, so no need to repeat it each time
    if (args$regressmt) {
        vars_to_regress <- c("mito_percentage")
    }
    if (args$ntgr == "none" | length(splitted_seurat_data) == 1){
        base::print(base::paste("Skipping datasets integration: either forced or only one identity is present."))
        if (args$gexnorm == "log"){
            base::print("Applying LogNormalize")
            scaled_norm_seurat_data <- Seurat::NormalizeData(splitted_seurat_data[[1]], verbose=FALSE)
            scaled_norm_seurat_data <- Seurat::FindVariableFeatures(
                scaled_norm_seurat_data,
                nfeatures=args$highvargex,
                verbose=FALSE
            )
            scaled_norm_seurat_data <- Seurat::ScaleData(
                scaled_norm_seurat_data,
                vars.to.regress=vars_to_regress,
                verbose=FALSE
            )
            SeuratObject::DefaultAssay(scaled_norm_seurat_data) <- "RNA"
        } else {
            base::print("Applying SCTransform")
            scaled_norm_seurat_data <- Seurat::SCTransform(
                splitted_seurat_data[[1]],
                assay="RNA",
                new.assay.name="SCT",
                variable.features.n=args$highvargex,
                vars.to.regress=vars_to_regress,
                verbose=FALSE                                                                      # too many stdout
            )
            SeuratObject::DefaultAssay(scaled_norm_seurat_data) <- "SCT"
        }
        base::rm(splitted_seurat_data)                                                             # remove unused data
        base::gc()
        return (scaled_norm_seurat_data)
    } else {
        base::print(base::paste("Running datasets integration"))
        if (args$gexnorm == "log"){
            base::print("Applying LogNormalize")
            for (i in 1:length(splitted_seurat_data)) {
                splitted_seurat_data[[i]] <- Seurat::NormalizeData(splitted_seurat_data[[i]], verbose=FALSE)
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
                k.weight=min(min(table(Idents(seurat_data))), 100),  # 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering
                verbose=FALSE
            )
            integrated_seurat_data <- Seurat::ScaleData(
                integrated_seurat_data,
                vars.to.regress=vars_to_regress,
                verbose=FALSE
            )
        } else {
            base::print("Applying SCTransform")
            for (i in 1:length(splitted_seurat_data)) {
                splitted_seurat_data[[i]] <- Seurat::SCTransform(
                    splitted_seurat_data[[i]],
                    assay="RNA",
                    new.assay.name="SCT",
                    variable.features.n=args$highvargex,
                    vars.to.regress=vars_to_regress,
                    verbose=FALSE                                                                      # too many stdout
                )
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
                k.weight=min(min(base::table(SeuratObject::Idents(seurat_data))), 100),  # 100 by default, but shouldn't be bigger than the min number of cells among all identities after filtering
                verbose=FALSE
            )
        }
        SeuratObject::DefaultAssay(integrated_seurat_data) <- "gex_integrated"
        base::rm(splitted_seurat_data, integration_features, integration_anchors)                  # remove unused data
        base::gc()
        return (integrated_seurat_data)
    }
}

gex_analyze <- function(seurat_data, args) {
    SeuratObject::Idents(seurat_data) <- "new.ident"                     # safety measure
    seurat_data <- gex_integrate(seurat_data, args)                      # sets "gex_integrated" as a default assay for integrated data, and either "RNA" or "SCT" for not integrated data
    base::print(
        base::paste(
            "Performing PCA reduction on", SeuratObject::DefaultAssay(seurat_data),
            "assay using 50 principal components"
        )
    )
    seurat_data <- Seurat::RunPCA(seurat_data, npcs=50, verbose=args$verbose)
    seurat_data <- Seurat::RunUMAP(
        seurat_data,
        reduction="pca",
        dims=args$gexndim,
        reduction.name="rnaumap",
        reduction.key="RNAUMAP_",
        verbose=args$verbose
    )
    return (seurat_data)
}