#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 1000 * 1024^2)  # 1GB should be good by default


suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(ggplot2))
suppressMessages(library(argparse))
suppressMessages(library(patchwork))


load_data <- function(location, mincells, cellnamedelim, cellnamefield) {
    raw_data <- Read10X(data.dir=location)
    return (
        CreateSeuratObject(
            counts=raw_data,
            min.cells=mincells,
            names.delim=cellnamedelim,
            names.field=cellnamefield
        )
    )
}


apply_qc_filters <- function(raw_seurat_data, minfeatures, maxmt, mtpattern, ribopattern) {
    qc_filtered_seurat_data <- PercentageFeatureSet(
        raw_seurat_data,
        pattern=mtpattern,
        col.name="percent_mt"
    )
    qc_filtered_seurat_data <- PercentageFeatureSet(
        qc_filtered_seurat_data,
        pattern=ribopattern,
        col.name="percent_ribo"
    )
    qc_filtered_seurat_data <- subset(
        qc_filtered_seurat_data,
        subset = nFeature_RNA >= minfeatures && percent_mt <= maxmt
    )
    return (qc_filtered_seurat_data)
}


export_vln_plot <- function(data, features, rootname, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(VlnPlot(data, features=features))
            dev.off()
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(VlnPlot(data, features=features))
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


export_pca_plot <- function(data, rootname, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(DimPlot(data, reduction="pca"))
            dev.off()
            pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
            print(DimPlot(data, reduction="pca"))
            dev.off()
            cat(paste("\nExport PCA plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export PCA plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_pca_heatmap <- function(data, rootname, dims=1:2, nfeatures=30, width=800, height=800, resolution=72){
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


export_pca_loadings_plot <- function(data, rootname, dims=1:2, nfeatures=30, width=800, height=800, resolution=72){
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

export_elbow_plot <- function(data, rootname, ndims=30, width=800, height=800, resolution=72){
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
    parser <- ArgumentParser(description='Run Seurat for Cellranger output')
    # Import data
    parser$add_argument("--mex",           help="Path to the folder with feature-barcode matrices in MEX format", type="character", required="True")
    parser$add_argument("--mincells",      help="Include features detected in at least this many cells. Default: 0", type="integer", default=0)
    parser$add_argument("--cellnamedelim", help="Delimiter to detect the initial type for each cell. Default: -", type="character", default="-")
    parser$add_argument("--cellnamefield", help="Field number to detect the initial type for each cell. Default: 2", type="integer", default=2)
    # Apply QC filters
    parser$add_argument("--minfeatures",   help="Include cells where at least this many features are detected. Default: 500", type="integer", default=500)
    parser$add_argument("--maxmt",         help="Maximum allowed mitochondrial contamination percentage. Default: 5", type="double", default=5)
    parser$add_argument("--highvarcount",  help="Number of higly variable features to detect. Default: 2000", type="integer", default=2000)
    parser$add_argument("--mtpattern",     help="Regex pattern to identify mitochondrial reads. Default: ^Mt-", type="character", default="^Mt-")
    parser$add_argument("--ribopattern",   help="Regex pattern to identify ribosomal reads. Default: ^Rp[sl]", type="character", default="^Rp[sl]")
    # Export results
    parser$add_argument("--output",        help='Output prefix. Default: ./seurat', type="character", default="./seurat")
    # Parallelization
    parser$add_argument("--threads",       help="Threads. Default: 1", type="integer", default=1)
    args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
    return (args)
}


args <- get_args()

print(paste("Setting parallelizations threads to", args$threads))
plan("multiprocess", workers=args$threads)
plan()

print(paste("Loading feature-barcode matrices from", args$mex))
raw_seurat_data <- load_data(args$mex, args$mincells, args$cellnamedelim, args$cellnamefield)
head(raw_seurat_data@meta.data, 5)

print("Exporting violin plots for raw data")
export_vln_plot(
    data=raw_seurat_data,
    features=c("nCount_RNA", "nFeature_RNA"),
    rootname=paste(args$output, "_raw_vln_plot", sep="")
)

print("Applying QC filters")
qc_filtered_seurat_data <- apply_qc_filters(
    raw_seurat_data,
    args$minfeatures,
    args$maxmt,
    args$mtpattern,
    args$ribopattern
)
head(qc_filtered_seurat_data@meta.data, 5)

print("Exporting violin plots for QC filtered data")
export_vln_plot(
    data=qc_filtered_seurat_data,
    features=c("nCount_RNA", "nFeature_RNA", "percent_mt", "percent_ribo"),
    rootname=paste(args$output, "_filtered_vln_plot", sep="")
)

print("Normalizing QC filtered data")
normalized_seurat_data <- NormalizeData(qc_filtered_seurat_data)
head(normalized_seurat_data@meta.data, 5)

print("Searching for the highly variable features")
with_highly_variable_features <- FindVariableFeatures(
    normalized_seurat_data,
    selection.method="vst",
    nfeatures=args$highvarcount
)
head(with_highly_variable_features)

print("Exporting highly variable feature plot")
export_variable_feature_plot(
    data=with_highly_variable_features,
    rootname=paste(args$output, "_highly_variable_feature_plot", sep="")
)

print("Scaling highly variable features")
scaled_data <- ScaleData(with_highly_variable_features)
head(scaled_data)

print("Performing linear dimensional reduction")
with_pca_data <- RunPCA(scaled_data)
head(with_pca_data)

print("Export Elbow plot to evaluate the dimensionality of the dataset")
export_elbow_plot(
    data=with_pca_data,
    rootname=paste(args$output, "_elbow_plot", sep="")
)

print("Export PCA, heatmap, and loadings plots for the first two principal components")
export_pca_plot(
    data=with_pca_data,
    rootname=paste(args$output, "_pca_plot", sep="")
)
export_pca_heatmap(
    data=with_pca_data,
    rootname=paste(args$output, "_pca_heatmap", sep="")
)
export_pca_loadings_plot(
    data=with_pca_data,
    rootname=paste(args$output, "_pca_loadings_plot", sep="")
)
