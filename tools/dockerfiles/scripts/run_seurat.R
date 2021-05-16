#!/usr/bin/env Rscript
options(warn=-1)
options("width"=200)
options(future.globals.maxSize = 4000 * 1024^2)  # 1GB should be good by default

suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(Seurat))
suppressMessages(library(future))
suppressMessages(library(ggplot2))
suppressMessages(library(argparse))
suppressMessages(library(patchwork))
suppressMessages(library(sctransform))


####################################################################################
# For details and treshold values look at
#    https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html
#    https://github.com/satijalab/seurat/issues/1679
####################################################################################


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


load_cell_cycle_data <- function (location) {
    cell_cycle_data <- read.table(
        location,
        sep=get_file_type(location),
        header=TRUE,
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
    return (cell_cycle_data)
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
            print(paste("Applying defaults - failed to load condition data from ", location))
            return (default_condition_data)
        }
    }
    print("Condition data is not provided. Applying defaults")
    return (default_condition_data)
}


load_seurat_data <- function(location, mincells, cell_identity_data, condition_data) {
    seurat_data <- CreateSeuratObject(
        counts=Read10X(data.dir=location),
        min.cells=mincells,
        names.delim="-",  # to get cell identity index from Cellranger aggr output
        names.field=2
    )
    seurat_data[["new.ident"]] <- cell_identity_data$library_id[Idents(seurat_data)]
    seurat_data[["condition"]] <- condition_data$condition[match(seurat_data$new.ident, condition_data$library_id)]
    Idents(seurat_data) <- "new.ident"
    return (seurat_data)
}


add_qc_metrics <- function(seurat_data, mitopattern) {
    # Number of genes detected per UMI: this metric with give us an idea of the
    # complexity of our dataset (more genes detected per UMI, more complex our data)
    # Mitochondrial percentage: this metric will give us a percentage of cell reads
    # originating from the mitochondrial genes
    seurat_data$log10_gene_per_log10_umi <- log10(seurat_data$nFeature_RNA) / log10(seurat_data$nCount_RNA)
    seurat_data$mito_percentage <- PercentageFeatureSet(seurat_data, pattern=mitopattern) 
    return (seurat_data)
}


apply_filters <- function(seurat_data, minfeatures, minumi, minnovelty, maxmt) {
    filtered_seurat_data <- subset(
        seurat_data,
        subset = (nFeature_RNA >= minfeatures) & 
                 (nCount_RNA >= minumi) &
                 (log10_gene_per_log10_umi >= minnovelty) &
                 (mito_percentage <= maxmt)
    )
    return (filtered_seurat_data)
}


explore_unwanted_variation <- function(seurat_data, cell_cycle_data, args) {
     temp_seurat_data <- NormalizeData(seurat_data, verbose=FALSE)
     temp_seurat_data <- CellCycleScoring(
         temp_seurat_data,
         s.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="s", "gene_id"]),
         g2m.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="g2/m", "gene_id"])
    )
    mito_quartiles <- quantile(temp_seurat_data@meta.data$mito_percentage, c(0.25, 0.5, 0.75))
    temp_seurat_data@meta.data$mito_factor <- cut(
        temp_seurat_data@meta.data$mito_percentage, 
        breaks=c(-Inf, mito_quartiles[1], mito_quartiles[2], mito_quartiles[3], Inf), 
        labels=c("Low", "Medium", "Medium high", "High")
    )
    temp_seurat_data <- FindVariableFeatures(temp_seurat_data, verbose=FALSE)
    temp_seurat_data <- ScaleData(temp_seurat_data, verbose=FALSE)
    temp_seurat_data <- RunPCA(temp_seurat_data, verbose=FALSE)
    temp_seurat_data <- RunUMAP(
        temp_seurat_data,
        reduction="pca",
        dims=1:args$ndim
    )
    export_dim_plot(
        data=temp_seurat_data,
        reduction="pca",
        split_by="Phase",
        group_by="Phase",
        rootname=paste(args$output, "_cell_cycle_effect_pca_plot", sep="")
    )
    export_dim_plot(
        data=temp_seurat_data,
        reduction="pca",
        split_by="mito_factor",
        group_by="mito_factor",
        rootname=paste(args$output, "_mitochondrial_gene_effect_pca_plot", sep="")
    )
    export_dim_plot(
        data=temp_seurat_data,
        reduction="umap",
        split_by="new.ident",
        rootname=paste(args$output, "_unintegrated_umap_plot_split_by_ident", sep="")
    )
}


integrate_seurat_data <- function(seurat_data, cell_cycle_data, regresscellcycle, regressmt, highvarcount) {
    splitted_seurat_data <- SplitObject(seurat_data, split.by="new.ident")
    for (i in 1:length(splitted_seurat_data)) {
        vars_to_regress <- NULL
        if (regressmt){
            vars_to_regress <- c("mito_percentage")
        }
        splitted_seurat_data[[i]] <- SCTransform(
            splitted_seurat_data[[i]],
            assay="RNA",
            new.assay.name="SCT",
            variable.features.n=highvarcount,
            vars.to.regress=vars_to_regress,
            verbose=FALSE
        )
        splitted_seurat_data[[i]] <- CellCycleScoring(
            splitted_seurat_data[[i]],
            s.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="s", "gene_id"]),
            g2m.features=as.vector(cell_cycle_data[tolower(cell_cycle_data$phase)=="g2/m", "gene_id"]),
            assay="SCT",
            verbose=FALSE
        )
        if (regresscellcycle){
            if (!is.null(vars_to_regress)) {
                vars_to_regress <- append(vars_to_regress, c("S.Score", "G2M.Score"))
            } else {
                vars_to_regress <- c("S.Score", "G2M.Score")
            }
            splitted_seurat_data[[i]] <- SCTransform(
                splitted_seurat_data[[i]],
                assay="RNA",
                new.assay.name="SCT",
                variable.features.n=highvarcount,
                vars.to.regress=vars_to_regress,
                verbose=FALSE
            )            
        }
    }
    integration_features <- SelectIntegrationFeatures(splitted_seurat_data, nfeatures=highvarcount)
    splitted_seurat_data <- PrepSCTIntegration(
        splitted_seurat_data, 
        anchor.features=integration_features,
        verbose=FALSE
    )
    integration_anchors <- FindIntegrationAnchors(
        splitted_seurat_data,
        normalization.method="SCT",
        anchor.features=integration_features,
        verbose=FALSE
    )
    integrated_seurat_data <- IntegrateData(
        integration_anchors, 
        new.assay.name="integrated",
        normalization.method="SCT",
        verbose=FALSE
    )
    DefaultAssay(integrated_seurat_data) <- "integrated"
    return (integrated_seurat_data)
}


get_conserved_markers <- function(cluster, seurat_data){
    conserved_markers <- FindConservedMarkers(
        seurat_data,
        ident.1=cluster,
        grouping.var="condition",
        verbose=FALSE
    ) %>% rownames_to_column(var="gene") %>% cbind(cluster_id=cluster, .)
    return (conserved_markers)
}


export_geom_bar_plot <- function(data, rootname, x_axis, color_by, x_label, y_label, legend_title, plot_title, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                geom_bar() +
                xlab(x_label) +
                ylab(y_label) +
                guides(fill=guide_legend(legend_title)) +
                ggtitle(plot_title)
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                print(
                    ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                    geom_bar() +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(fill=guide_legend(legend_title)) +
                    ggtitle(plot_title)
                )
                dev.off()
            }
            cat(paste("\nExport bar plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export bar plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        }
    )
}

export_log10_geom_density_plot <- function(data, rootname, x_axis, color_by, x_intercept, x_label, y_label, legend_title, plot_title, alpha=0.2, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                geom_density(alpha=alpha) +
                scale_x_log10() + 
                xlab(x_label) +
                ylab(y_label) +
                guides(fill=guide_legend(legend_title)) +
                geom_vline(xintercept=x_intercept, color="red") +
                ggtitle(plot_title)
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                print(
                    ggplot(data, aes_string(x=x_axis, fill=color_by)) +
                    geom_density(alpha=alpha) +
                    scale_x_log10() + 
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(fill=guide_legend(legend_title)) +
                    geom_vline(xintercept=x_intercept, color="red") +
                    ggtitle(plot_title)
                )
                dev.off()
            }
            cat(paste("\nExport log10 density plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export log10 density plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        }
    )
}


export_log10_geom_point_plot <- function(data, rootname, x_axis, y_axis, color_by, facet_by, x_intercept, y_intercept, x_label, y_label, legend_title, plot_title, alpha=0.2, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                ggplot(data, aes_string(x=x_axis, y=y_axis, color=color_by)) +
                geom_point(alpha=alpha) +
                scale_colour_gradient(low="gray90", high="black") +
                stat_smooth(method=lm) +
                scale_x_log10() + 
  	            scale_y_log10() + 
                geom_vline(xintercept=x_intercept, color="red") +
  	            geom_hline(yintercept=y_intercept, color="red") +
                facet_wrap(as.formula(paste("~", facet_by))) +
                xlab(x_label) +
                ylab(y_label) +
                guides(color=guide_legend(legend_title)) +
                ggtitle(plot_title)
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                print(
                    ggplot(data, aes_string(x=x_axis, y=y_axis, color=color_by)) +
                    geom_point(alpha=alpha) +
                    scale_colour_gradient(low="gray90", high="black") +
                    stat_smooth(method=lm) +
                    scale_x_log10() + 
                    scale_y_log10() + 
                    geom_vline(xintercept=x_intercept, color="red") +
                    geom_hline(yintercept=y_intercept, color="red") +
                    facet_wrap(as.formula(paste("~", facet_by))) +
                    xlab(x_label) +
                    ylab(y_label) +
                    guides(color=guide_legend(legend_title)) +
                    ggtitle(plot_title)
                )
                dev.off()
            }
            cat(paste("\nExport log10 point plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export log10 point plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        }
    )
}


export_vln_plot <- function(data, features, rootname, group_by=NULL, pt_size=NULL, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                VlnPlot(
                    data,
                    features=features,
                    pt.size = pt_size,
                    group.by=group_by
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                print(
                    VlnPlot(
                        data,
                        features=features,
                        pt.size = pt_size,
                        group.by=group_by
                    )
                )
                dev.off()
            }
            cat(paste("\nExport violin plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export violin plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_feature_plot <- function(data, features, rootname, reduction, split_by=NULL, label=FALSE, order=FALSE, min_cutoff=NA, max_cutoff=NA, pt_size=NULL, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                FeaturePlot(
                    data,
                    features=features,
                    pt.size=pt_size,
                    order=order,
                    min.cutoff=min_cutoff,
                    max.cutoff=max_cutoff,
                    reduction=reduction,
                    split.by=split_by,
                    label=label
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                print(
                    FeaturePlot(
                        data,
                        features=features,
                        pt.size=pt_size,
                        order=order,
                        min.cutoff=min_cutoff,
                        max.cutoff=max_cutoff,
                        reduction=reduction,
                        split.by=split_by,
                        label=label
                    )
                )
                dev.off()
            }
            cat(paste("\nExport feature plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export feature plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_variable_feature_plot <- function(data, rootname, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(VariableFeaturePlot(data))
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                print(VariableFeaturePlot(data))
                dev.off()
            }
            cat(paste("\nExport variable feature plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export variable feature plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_dim_plot <- function(data, rootname, reduction, split_by=NULL, group_by=NULL, label=FALSE, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                DimPlot(
                    data,
                    reduction=reduction,
                    split.by=split_by,
                    group.by=group_by,
                    label=label
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                print(
                    DimPlot(
                        data,
                        reduction=reduction,
                        split.by=split_by,
                        group.by=group_by,
                        label=label
                    )
                )
                dev.off()
            }
            cat(paste("\nExport Dim plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export Dim plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_pca_heatmap <- function(data, rootname, dims=NULL, nfeatures=30, cells=500, pdf=FALSE, width=800, height=800, resolution=72){
    tryCatch(
        expr = {
            png(filename=paste(rootname, ".png", sep=""), width=width, height=height, res=resolution)
            print(
                DimHeatmap(
                    data,
                    dims=dims,
                    nfeatures=nfeatures,
                    reduction="pca",
                    cells=cells,
                    balanced=TRUE
                )
            )
            dev.off()
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                print(
                    DimHeatmap(
                        data,
                        dims=dims,
                        nfeatures=nfeatures,
                        reduction="pca",
                        cells=cells,
                        balanced=TRUE
                    )
                )
                dev.off()
            }
            cat(paste("\nExport PCA heatmap to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export PCA heatmap to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_pca_loadings_plot <- function(data, rootname, dims=NULL, nfeatures=30, pdf=FALSE, width=800, height=800, resolution=72){
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
            if (pdf) {
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
            }
            cat(paste("\nExport PCA loadings plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export PCA loadings plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_elbow_plot <- function(data, rootname, ndims=NULL, pdf=FALSE, width=800, height=800, resolution=72){
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
            if (pdf) {
                pdf(file=paste(rootname, ".pdf", sep=""), width=round(width/resolution), height=round(height/resolution))
                print(
                    ElbowPlot(
                        data,
                        ndims=ndims,
                        reduction="pca"
                    )
                )
                dev.off()
            }
            cat(paste("\nExport Elbow plot to ", rootname, ".(png/pdf)", "\n", sep=""))
        },
        error = function(e){
            dev.off()
            cat(paste("\nFailed to export Elbow plot to ", rootname, ".(png/pdf)", "\n",  sep=""))
        }
    )
}


export_all_clustering_plots <- function(seurat_data, suffix, args) {
    for (i in 1:length(args$resolution)) {
        current_resolution <- args$resolution[i]
        export_dim_plot(
            data=seurat_data,
            reduction="umap",
            split_by="condition",
            group_by=paste("integrated_snn_res", current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "umap_plot_split_by_condition_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        export_dim_plot(
            data=seurat_data,
            reduction="umap",
            split_by="Phase",
            group_by=paste("integrated_snn_res", current_resolution, sep="."),
            label=TRUE,
            rootname=paste(args$output, suffix, "umap_plot_split_by_phase_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        Idents(seurat_data) <- paste("integrated_snn_res", current_resolution, sep=".")
        export_feature_plot(
            data=seurat_data,
            features=c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "mito_percentage", "log10_gene_per_log10_umi"),
            reduction="umap",
            label=TRUE,
            order=TRUE,
            min_cutoff="q10",
            rootname=paste(args$output, suffix, "umap_feature_plot_res", current_resolution, sep="_"),
            pdf=args$pdf
        )
        Idents(seurat_data) <- "new.ident"
    }
}


export_all_dimensionality_plots <- function(seurat_data, suffix, args) {
    export_elbow_plot(
        data=seurat_data,
        ndims=50,
        rootname=paste(args$output, suffix, "elbow_plot", sep="_"),
        pdf=args$pdf
    )
    print("Export PCA plot for clustered integrated data")
    export_dim_plot(
        data=seurat_data,
        reduction="pca",
        rootname=paste(args$output, suffix, "pca_plot", sep="_"),
        pdf=args$pdf
    )
    print("Export heatmaps and loadings plots for all 50 principal components")
    export_pca_heatmap(
        data=seurat_data,
        dims=1:50,
        rootname=paste(args$output, suffix, "pca_heatmap", sep="_"),
        pdf=args$pdf
    )
    export_pca_loadings_plot(
        data=seurat_data,
        dims=1:50,
        rootname=paste(args$output, suffix, "pca_loadings_plot", sep="_"),
        pdf=args$pdf
    )
    export_dim_plot(
        data=seurat_data,
        reduction="umap",
        split_by="new.ident",
        rootname=paste(args$output, suffix, "umap_plot_split_by_ident", sep="_"),
        pdf=args$pdf
    )
}


export_all_qc_plots <- function(seurat_data, suffix, args){
    export_geom_bar_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "cell_count_plot", sep="_"),
        x_axis="new.ident",
        color_by="new.ident",
        x_label="Identity",
        y_label="Cells",
        legend_title="Identity",
        plot_title=paste("Number of cells per dataset (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_log10_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "umi_log10_density_plot", sep="_"),
        x_axis="nCount_RNA",
        color_by="new.ident",
        x_intercept=args$minumi,
        x_label="UMI's per cell",
        y_label="Density (log10)",
        legend_title="Identity",
        plot_title=paste("UMI density per cell (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_log10_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gene_log10_density_plot", sep="_"),
        x_axis="nFeature_RNA",
        color_by="new.ident",
        x_intercept=args$minfeatures,
        x_label="Genes per cell",
        y_label="Density (log10)",
        legend_title="Identity",
        plot_title=paste("Gene density per cell (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_log10_geom_point_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "gene_umi_log10_correlation_plot", sep="_"),
        x_axis="nCount_RNA",
        y_axis="nFeature_RNA",
        color_by="mito_percentage",
        facet_by="new.ident",
        x_intercept=args$minumi,
        y_intercept=args$minfeatures,
        x_label="UMI's (log10)",
        y_label="Genes (log10)",
        legend_title="Mitochondrial %",
        plot_title=paste("Gene vs UMI correlation (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_log10_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "mitochondrial_gene_percentage_log10_density_plot", sep="_"),
        x_axis="mito_percentage",
        color_by="new.ident",
        x_intercept=args$maxmt,
        x_label="Mitochondrial gene percentage per cell",
        y_label="Density (log10)",
        legend_title="Identity",
        plot_title=paste("Mitochondrial gene percentage density per cell (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_log10_geom_density_plot(
        data=seurat_data@meta.data,
        rootname=paste(args$output, suffix, "novelty_score_log10_density_plot", sep="_"),
        x_axis="log10_gene_per_log10_umi",
        color_by="new.ident",
        x_intercept=args$minnovelty,
        x_label="log10 Gene / log10 UMI per cell",
        y_label="Density (log10)",
        legend_title="Identity",
        plot_title=paste("Novelty score (log10Gene/log10UMI) density per cell (", suffix, ")", sep=""),
        pdf=args$pdf
    )
    export_vln_plot(
        data=seurat_data,
        features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi"),
        rootname=paste(args$output, suffix, "vln_plot", sep="_"),
        pt_size=0,
        pdf=args$pdf
    )
    export_vln_plot(
        data=seurat_data,
        features=c("nCount_RNA", "nFeature_RNA", "mito_percentage", "log10_gene_per_log10_umi"),
        rootname=paste(args$output, suffix, "vln_plot_grouped_by_condition", sep="_"),
        group_by="condition",
        pt_size=0,
        pdf=args$pdf
    )
}


get_args <- function(){
    parser <- ArgumentParser(description='Runs Seurat for comparative scRNA-seq analysis of across experimental conditions')
    # Import data from Cellranger Aggregate results
    parser$add_argument("--mex",           help="Path to the folder with not normalized not filtered aggregated feature-barcode matrices in MEX format", type="character", required="True")
    parser$add_argument("--identity",      help="Path to the aggregation CSV file to set the initial cell identity classes", type="character", required="True")
    parser$add_argument("--cellcycle",     help="Path to the TSV/CSV file with cell cycle data. First column - 'phase', second column 'gene_id'", type="character", required="True")
    parser$add_argument("--condition",     help="Path to the TSV/CSV file to define datasets conditions for grouping. First column - 'library_id' with values from the --identity file, second column 'condition'. Default: each dataset is assigned to its own biological condition", type="character")    
    # Apply QC filters
    parser$add_argument("--mincells",      help="Include features detected in at least this many cells (applied to thoughout all datasets together). Default: 10", type="integer", default=10)
    parser$add_argument("--minfeatures",   help="Include cells where at least this many features are detected. Default: 250", type="integer", default=250)
    parser$add_argument("--minumi",        help="Include cells where at least this many UMI are detected. Default: 500", type="integer", default=500)
    parser$add_argument("--minnovelty",    help="Include cells with the novelty score not lower that this value (calculated as log10(genes)/log10(UMIs)). Default: 0.8", type="double", default=0.8)
    parser$add_argument("--maxmt",         help="Include cells with the mitochondrial contamination percentage not bigger that this value. Default: 5", type="double", default=5)
    parser$add_argument("--mitopattern",   help="Regex pattern to identify mitochondrial reads. Default: ^Mt-", type="character", default="^Mt-")
    # Integration and clustering parameters
    parser$add_argument("--regresscellcycle", help="Regress cell cycle as a confounding sources of variation. Default: false", action="store_true")
    parser$add_argument("--regressmt",     help="Regress mitochondrial gene expression as a confounding sources of variation. Default: false", action="store_true")
    parser$add_argument("--highvarcount",  help="Number of higly variable features to detect. Default: 3000", type="integer", default=3000)
    parser$add_argument("--ndim",          help="Number of principal components to use in clustering (1:50). Use Elbow plot to adjust this parameter. Default: 10", type="integer", default=10)
    parser$add_argument("--resolution",    help="Clustering resolution. Can be set as array. Default: 0.4 0.6 0.8 1.0 1.4", type="double", default=c(0.4, 0.6, 0.8, 1.0, 1.4), nargs="*")
    # Export results
    parser$add_argument("--pdf",           help="Export plots in PDF. Default: false", action="store_true")
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
print(paste("Loading cell cycle data from", args$cellcycle))
cell_cycle_data <- load_cell_cycle_data(args$cellcycle)
print("Trying to load condition data")
condition_data <- load_condition_data(args$condition, cell_identity_data)
print(paste("Loading feature-barcode matrices from", args$mex))
seurat_data <- load_seurat_data(args$mex, args$mincells, cell_identity_data, condition_data)
print("Adding QC metrics to not filtered seurat data")
seurat_data <- add_qc_metrics(seurat_data, args$mitopattern)
export_all_qc_plots(seurat_data, "raw", args)
print("Applying QC filters to all datasets at once")
seurat_data <- apply_filters(seurat_data, args$minfeatures, args$minumi, args$minnovelty, args$maxmt)
export_all_qc_plots(seurat_data, "filtered", args)
print("Evaluating effects of cell cycle and mitochodrial gene expression")
explore_unwanted_variation(seurat_data, cell_cycle_data, args)
print("Running dataset integration")
seurat_data <- integrate_seurat_data(seurat_data, cell_cycle_data, args$regresscellcycle, args$regressmt, args$highvarcount)
head(seurat_data@meta.data)
print("Performing PCA reduction of integrated data. Use all 50 principal components")
seurat_data <- RunPCA(seurat_data, npcs=50, verbose=FALSE)
print(paste("Performing UMAP reduction of integrated data using", args$ndim, "principal components"))
seurat_data <- RunUMAP(seurat_data, reduction="pca", dims=1:args$ndim)
export_all_dimensionality_plots(seurat_data, "integrated", args)
print(paste("Clustering integrated data using", args$ndim, "principal components"))
seurat_data <- FindNeighbors(seurat_data, reduction="pca", dims=1:args$ndim)
seurat_data <- FindClusters(seurat_data, resolution=args$resolution)
export_all_clustering_plots(seurat_data, "integrated_clustered", args)




print("Identifying conserved markers for all clusters irrespective of condition")
DefaultAssay(seurat_data) <- "RNA"
conserved_markers <- map_dfr(
    sort(unique(seurat_data@meta.data$seurat_clusters)),  # need to get all clusters
    seurat_data,
    get_conserved_markers
)
head(conserved_markers)