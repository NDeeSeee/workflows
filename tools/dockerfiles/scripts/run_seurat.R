#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 3000 * 1024^2)  # 1GB should be good by default


suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(ggplot2))
suppressMessages(library(argparse))
suppressMessages(library(patchwork))


set_threads <- function (threads) {
    plan("multiprocess", workers=threads)
    plan()
}


get_file_type <- function (filename) {
    ext = tools::file_ext(filename)
    separator = "\t"
    if (ext == "csv"){
        separator = ","
    }
    return (separator)
}


load_cell_identity_data <- function (location) {
    cell_identity_data <- read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    return (cell_identity_data)
}


load_condition_data <- function(location, cell_identity_data) {
    default_condition_data <- data.frame(
        library_id=cell_identity_data$library_id,
        condition=rownames(cell_identity_data),
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
        if (all(is.element(cell_identity_data$library_id, condition_data$library_id))){
            print(paste("Condition data is successfully loaded from ", location))
            return (condition_data)
        } else {
            print("Condition data is not well formatted. Applying defaults")
            return (default_condition_data)
        }
    }
    print("Condition data is not provided. Applying defaults")
    return (default_condition_data)
}


load_seurat_data <- function(location, cell_identity_data, condition_data) {
    seurat_data <- CreateSeuratObject(
        counts=Read10X(data.dir=location),
        names.delim="-",  # to get cell identity index from Cellranger aggr output
        names.field=2
    )
    seurat_data[["new.ident"]] <- cell_identity_data$library_id[Idents(seurat_data)]
    seurat_data[["condition"]] <- condition_data$condition[match(seurat_data$new.ident, condition_data$library_id)]
    Idents(seurat_data) <- "new.ident"
    return (seurat_data)
}


apply_filters <- function(raw_seurat_data, minfeatures, maxmt, mtpattern) {
    filtered_seurat_data <- PercentageFeatureSet(
        raw_seurat_data,
        pattern=mtpattern,
        col.name="percent_mt"
    )
    filtered_seurat_data <- subset(
        filtered_seurat_data,
        subset = nFeature_RNA >= minfeatures && percent_mt <= maxmt
    )
    return (filtered_seurat_data)
}


integrate_seurat_data <- function(seurat_data, highvarcount) {
    print("Splitting filtered data by condition")
    splitted_seurat_data <- SplitObject(seurat_data, split.by="condition")
    print("Normalizing and identifying variable features for each condition")
    for (i in 1:length(splitted_seurat_data)) {
        condition <- splitted_seurat_data[[i]]@meta.data$condition[1]
        print(paste("Processing condition", condition))
        splitted_seurat_data[[i]] <- NormalizeData(
            splitted_seurat_data[[i]],
            verbose=FALSE
        )
        splitted_seurat_data[[i]] <- FindVariableFeatures(
            splitted_seurat_data[[i]],
            selection.method="vst",
            nfeatures=highvarcount,
            verbose=FALSE
        )
        export_variable_feature_plot(
            data=splitted_seurat_data[[i]],
            rootname=paste(args$output, condition, "highly_variable_feature_plot", sep="_")
        )
    }
    print("Selecting features that are repeatedly variable accross datasets for integration")
    integration_features <- SelectIntegrationFeatures(splitted_seurat_data)
    print("Searching for integration anchors")
    integration_anchors <- FindIntegrationAnchors(
        splitted_seurat_data,
        anchor.features=integration_features
    )
    print("Creating integrated data assay")
    integrated_seurat_data <- IntegrateData(
        integration_anchors, 
        new.assay.name="integrated",
        verbose=FALSE
    )
    DefaultAssay(integrated_seurat_data) <- "integrated"
    return (integrated_seurat_data)
}


export_vln_plot <- function(data, features, rootname, group_by=NULL, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                VlnPlot(
                    data,
                    features=features,
                    group.by=group_by
                )
            )
            dev.off()
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(
                VlnPlot(
                    data,
                    features=features,
                    group.by=group_by
                )
            )
            dev.off()
            cat(paste("\nExport violin plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export violin plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_variable_feature_plot <- function(data, rootname, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(VariableFeaturePlot(data))
            dev.off()
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(VariableFeaturePlot(data))
            dev.off()
            cat(paste("\nExport variable feature plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export variable feature plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_dim_plot <- function(data, rootname, reduction, split_by=NULL, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                DimPlot(
                    data,
                    reduction=reduction,
                    split.by=split_by
                )
            )
            dev.off()
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(
                DimPlot(
                    data,
                    reduction=reduction,
                    split.by=split_by
                )
            )
            dev.off()
            cat(paste("\nExport Dim plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export Dim plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_pca_heatmap <- function(data, rootname, dims=NULL, nfeatures=30, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                DimHeatmap(
                    data,
                    dims=dims,
                    nfeatures=nfeatures,
                    reduction="pca",
                    balanced=TRUE
                )
            )
            dev.off()
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(
                DimHeatmap(
                    data,
                    dims=dims,
                    nfeatures=nfeatures,
                    reduction="pca",
                    balanced=TRUE
                )
            )
            dev.off()
            cat(paste("\nExport PCA heatmap to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export PCA heatmap to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_pca_loadings_plot <- function(data, rootname, dims=NULL, nfeatures=30, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                VizDimLoadings(
                    data,
                    dims=dims,
                    nfeatures=nfeatures,
                    reduction="pca"
                )
            )
            dev.off()
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(
                VizDimLoadings(
                    data,
                    dims=dims,
                    nfeatures=nfeatures,
                    reduction="pca"
                )
            )
            dev.off()
            cat(paste("\nExport PCA loadings plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export PCA loadings plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}

export_elbow_plot <- function(data, rootname, ndims=NULL, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                ElbowPlot(
                    data,
                    ndims=ndims,
                    reduction="pca"
                )
            )
            dev.off()
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(
                ElbowPlot(
                    data,
                    ndims=ndims,
                    reduction="pca"
                )
            )
            dev.off()
            cat(paste("\nExport Elbow plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export Elbow plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


get_args <- function(){
    parser <- ArgumentParser(description='Runs Seurat for comparative scRNA-seq analysis of across experimental conditions')
    # Import data from Cellranger Aggregate results
    parser$add_argument("--mex",           help="Path to the folder with aggregated feature-barcode matrices in MEX format", type="character", required="True")
    parser$add_argument("--identity",      help="Path to the aggregation CSV file to set the initial cell identity classes", type="character", required="True")
    parser$add_argument("--condition",     help="Path to the TSV/CSV file to define datasets conditions for grouping. First column - 'library_id' with values from the --identity file, second column 'condition'. Default: each dataset is assigned to its own biological condition", type="character")
    # Apply QC filters
    parser$add_argument("--minfeatures",   help="Include cells where at least this many features are detected (per sample). Default: 500", type="integer", default=500)
    parser$add_argument("--maxmt",         help="Maximum allowed mitochondrial contamination percentage (per sample). Default: 5", type="double", default=5)
    parser$add_argument("--mtpattern",     help="Regex pattern to identify mitochondrial reads. Default: ^Mt-", type="character", default="^Mt-")
    # Integration and clustering parameters
    parser$add_argument("--highvarcount",  help="Number of higly variable features to detect. Default: 2000", type="integer", default=2000)
    parser$add_argument("--ndim",          help="Number of principal components to use in clustering (1:50). Use Elbow plot to adjust this parameter. Default: 10", type="integer", default=10)
    parser$add_argument("--resolution",    help="Clustering resolution. Default: 0.8", type="double", default=0.8)
    # Export results
    parser$add_argument("--output",        help="Output prefix. Default: ./seurat", type="character", default="./seurat")
    # Performance parameters
    parser$add_argument("--threads",       help="Threads. Default: 1", type="integer", default=1)
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}


args <- get_args()

print(paste("Setting parallelizations threads to", args$threads))
set_threads(args$threads)

print(paste("Loading cell identity data from Cellranger aggregation metadata file", args$identity))
cell_identity_data <- load_cell_identity_data (args$identity)
head(cell_identity_data)

print("Loading condition data or setting it to default values")
condition_data <- load_condition_data(args$condition, cell_identity_data)
head(condition_data)

print(paste("Loading seurat data from feature-barcode matrices in", args$mex))
seurat_data <- load_seurat_data(args$mex, cell_identity_data, condition_data)
head(seurat_data@meta.data)

print("Exporting violin plots for not filtered seurat data")
export_vln_plot(
    data=seurat_data,
    features=c("nCount_RNA", "nFeature_RNA"),
    rootname=paste(args$output, "_raw_vln_plot", sep="")
)
export_vln_plot(
    data=seurat_data,
    features=c("nCount_RNA", "nFeature_RNA"),
    group_by="condition",
    rootname=paste(args$output, "_raw_grouped_by_condition_vln_plot", sep="")
)

print("Applying QC filters to all datasets at once")  # These filters are applied to each cell individually, so we can apply them for all datasets at once
seurat_data <- apply_filters(
    seurat_data,
    args$minfeatures,
    args$maxmt,
    args$mtpattern
)
head(seurat_data@meta.data)

print("Exporting violin plots for QC filtered data")
export_vln_plot(
    data=seurat_data,
    features=c("nCount_RNA", "nFeature_RNA", "percent_mt"),
    rootname=paste(args$output, "_filtered_vln_plot", sep="")
)
export_vln_plot(
    data=seurat_data,
    features=c("nCount_RNA", "nFeature_RNA", "percent_mt"),
    group_by="condition",
    rootname=paste(args$output, "_filtered_grouped_by_condition_vln_plot", sep="")
)

print("Running dataset integration splitting by condition")
integrated_seurat_data <- integrate_seurat_data(
    seurat_data,
    args$highvarcount
)
head(integrated_seurat_data@meta.data)

print("Scaling integrated data")
integrated_seurat_data <- ScaleData(integrated_seurat_data, verbose=FALSE)
head(integrated_seurat_data)

print("Performing PCA reduction of integrated data. Use all 50 principal components")
integrated_seurat_data <- RunPCA(
    integrated_seurat_data,
    npcs=50,
    verbose=FALSE
)
head(integrated_seurat_data)

print("Export Elbow plot to evaluate the dimensionality of the dataset. Use all 50 principal components")
export_elbow_plot(
    data=integrated_seurat_data,
    ndims=50,
    rootname=paste(args$output, "_elbow_plot", sep="")
)

print("Export PCA plot for clustered integrated data")
export_dim_plot(
    data=integrated_seurat_data,
    reduction="pca",
    split_by="condition",
    rootname=paste(args$output, "_pca_plot", sep="")
)
print(paste("Export heatmaps and loadings plots for", args$ndim, "principal components"))
export_pca_heatmap(
    data=integrated_seurat_data,
    dims=1:args$ndim,
    rootname=paste(args$output, "_pca_heatmap", sep="")
)
export_pca_loadings_plot(
    data=integrated_seurat_data,
    dims=1:args$ndim,
    rootname=paste(args$output, "_pca_loadings_plot", sep="")
)

print(paste("Performing UMAP reduction of integrated data using", args$ndim, "principal components"))
integrated_seurat_data <- RunUMAP(
    integrated_seurat_data,
    reduction="pca",
    dims=1:args$ndim
)
head(integrated_seurat_data)

print(paste("Clustering integrated data using", args$ndim, "principal components"))
integrated_seurat_data <- FindNeighbors(
    integrated_seurat_data,
    reduction="pca",
    dims=1:args$ndim
)
integrated_seurat_data <- FindClusters(
    integrated_seurat_data,
    resolution=args$resolution
)

print("Export UMAP plots for clustered integrated data")
export_dim_plot(
    data=integrated_seurat_data,
    reduction="umap",
    split_by="condition",
    rootname=paste(args$output, "_integrated_umap_plot", sep="")
)