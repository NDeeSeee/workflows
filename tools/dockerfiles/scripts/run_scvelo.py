import os
import sys
import argparse
import anndata
import scvelo
import pandas
import numpy
import matplotlib


def configure_scvelo(args):
    scvelo.set_figure_params(
        style="scvelo",
        dpi_save=args.dpi,
        transparent=False
    )
    scvelo.settings.verbosity = 3
    scvelo.settings.plot_prefix = args.output
    scvelo.settings.figdir = ""
    scvelo.settings.autoshow = False


def get_velocity_data(args):
    print(f"""velocity - {args.h5ad}\nbarcodes - {args.cells}\nUMAP     - {args.umap}\nPCA      - {args.pca}\nclusters - {args.clusters}""")
    raw_velocity_data = anndata.read_h5ad(args.h5ad)
    cells_data = pandas.read_csv(args.cells, header=None, sep="\t")
    filtered_velocity_data = raw_velocity_data[cells_data.iloc[:, 0]].copy()
    clusters_data = pandas.read_csv(args.clusters, header=None, sep="\t")
    filtered_velocity_data.obs["clusters"] = clusters_data.values
    try:
        # https://github.com/theislab/scvelo/issues/37#issuecomment-499101214
        # https://github.com/theislab/scvelo/discussions/722
        pca_data = pandas.read_csv(args.pca, header=None, sep="\t")
        umap_data = pandas.read_csv(args.umap, header=None, sep="\t")
        filtered_velocity_data.obsm["X_pca"] = pca_data.values
        filtered_velocity_data.obsm["X_umap"] = umap_data.values
    except Exception:
        print("Failed to load custom PCA/UMAP coordinates. Will be calculated internally")
    return filtered_velocity_data


def estimate_velocity(velocity_data, args):
    scvelo.pp.filter_and_normalize(velocity_data)
    if ("X_umap" not in velocity_data.obsm) or ("X_pca" not in velocity_data.obsm):
        print(f"""Computing moments with calculated PCA and UMAP using {args.ndim} PCs and {args.nneighbors} nearest neighbors""")
        scvelo.pp.moments(velocity_data, n_pcs=args.ndim, n_neighbors=args.nneighbors)
        scvelo.tl.umap(velocity_data)
    else:
        print(f"""Computing moments using custom PCA and UMAP coordinates and {args.nneighbors} nearest neighbors""")
        scvelo.pp.moments(velocity_data, n_pcs=None, n_neighbors=args.nneighbors)
    scvelo.tl.recover_dynamics(velocity_data, n_jobs=args.threads)
    scvelo.tl.velocity(velocity_data, mode="dynamical")
    scvelo.tl.velocity_graph(velocity_data, n_jobs=args.threads)
    scvelo.tl.velocity_confidence(velocity_data)
    scvelo.tl.rank_velocity_genes(velocity_data, groupby="clusters", min_corr=0.3, n_genes=args.top)
    scvelo.tl.velocity_pseudotime(velocity_data)
    scvelo.tl.rank_dynamical_genes(velocity_data, groupby="clusters", n_genes=args.top)
    scvelo.tl.latent_time(velocity_data)
    velocity_data.uns["neighbors"]["distances"] = velocity_data.obsp["distances"]
    velocity_data.uns["neighbors"]["connectivities"] = velocity_data.obsp["connectivities"]
    scvelo.tl.paga(
        velocity_data,
        groups="clusters",
        use_time_prior="latent_time"    # "velocity_pseudotime"
    )


def export_text_data(velocity_data, args):
    velocity_genes = scvelo.get_df(velocity_data, "rank_velocity_genes/names")
    velocity_genes.to_csv(args.output+"vel_genes.tsv", sep="\t")
    dynamical_genes = scvelo.get_df(velocity_data, "rank_dynamical_genes/names")
    dynamical_genes.to_csv(args.output+"dyn_genes.tsv", sep="\t")


def export_velocity_plot(velocity_data, args):
    scvelo.pl.proportions(
        velocity_data,
        groupby="clusters",
        show=False,              # need to specify explicitely for scvelo.pl.proportions
        save="proportions.png"
    )
    scvelo.pl.scatter(
        velocity_data,
        c=["velocity_length", "velocity_confidence"],
        cmap="coolwarm",
        perc=[5, 95],
        save="confidence.png"
    )
    scvelo.pl.velocity_embedding_stream(
        velocity_data,
        basis="umap",
        color="clusters",
        save="stream.png"
    )
    scvelo.pl.velocity_embedding_grid(
        velocity_data,
        basis="umap",
        color="clusters",
        save="grid.png"
    )
    scvelo.pl.scatter(
        velocity_data,
        color="velocity_pseudotime",
        alpha=0.25,
        cmap="gnuplot",
        save="pseudotime.png"
    )
    scvelo.pl.scatter(
        velocity_data,
        color="latent_time",
        color_map="gnuplot",
        alpha=0.25,
        save="latenttime.png"
    )
    scvelo.pl.paga(
        velocity_data,
        basis="umap",
        alpha=0.1,
        min_edge_width=2,
        node_size_scale=1.5,
        size=50,
        legend_loc="on data",
        legend_fontsize=6,
        fontsize=6,
        save="paga.png"
    )
    top_genes = velocity_data.var["fit_likelihood"].sort_values(ascending=False).index[:args.top]
    top_genes = top_genes.tolist()
    scvelo.pl.heatmap(
        velocity_data,
        var_names=top_genes,
        sortby="latent_time",
        col_color="clusters",
        n_convolve=300,
        font_scale=0.2,
        yticklabels=True,                                        # to display all genes
        save=f"""top_{args.top}_driver_genes_heatmap.png"""
    )
    # scvelo.pl.scatter(
    #     velocity_data,
    #     x="latent_time",
    #     y=top_genes,
    #     color="clusters",
    #     frameon=False,
    #     ncols=5,
    #     legend_loc="best",
    #     save=f"""top_{args.top}_driver_genes_scatter.png"""
    # )
    scvelo.pl.velocity(
        velocity_data,
        var_names=top_genes,
        color="clusters",
        save="driver_genes_velocity.png"
    )
    # if user provided genes of interest
    if args.genes and len(args.genes) > 0:
        scvelo.pl.velocity(
            velocity_data,
            var_names=args.genes,
            save="selected_genes_velocity.png"
        )


def normalize_args(args, skip_list=[]):
    normalized_args = {}
    for key,value in args.__dict__.items():
        if key not in skip_list:
            normalized_args[key] = value if not value or os.path.isabs(value) else os.path.normpath(os.path.join(os.getcwd(), value))
        else:
            normalized_args[key]=value
    return argparse.Namespace (**normalized_args)


def arg_parser():
    general_parser = argparse.ArgumentParser()
    general_parser.add_argument(
        "--h5ad",
        help="Path to the h5ad file",
        type=str,
        required=True
    )
    general_parser.add_argument(
        "--cells",
        help="Path to the cells data TSV file",
        type=str,
        required=True
    )
    general_parser.add_argument(
        "--umap",
        help="Path to the custom UMAP data TSV file. \
              Ignored if --pca is not provided",
        type=str
    )
    general_parser.add_argument(
        "--pca",
        help="Path to the custom PCA data TSV file.\
              Ignored if --umap is not provided",
        type=str
    )
    general_parser.add_argument(
        "--ndim",
        help="Number of principal components to use. \
              Ignored when both --pca and --umap are provided",
        type=int,
        default=30
    )
    general_parser.add_argument(
        "--nneighbors",
        help="Number of nearest neighbors to use when computing moments for velocity estimation",
        type=int,
        default=30
    )
    general_parser.add_argument(
        "--genes",
        help="List of genes of interest",
        type=str,
        nargs="*"
    )
    general_parser.add_argument(
        "--clusters",
        help="Path to the clusters data TSV file",
        type=str,
        required=True
    )
    general_parser.add_argument(
        "--dpi",
        help="Resolution for exported plots",
        type=int,
        default=300
    )
    general_parser.add_argument(
        "--threads",
        help="Number of jobs to run in parallel",
        type=int,
        default=1
    )
    general_parser.add_argument(
        "--top",
        help="Top N genes to display on plots and return in tables",
        type=int,
        default=100
    )
    general_parser.add_argument(
        "--output",
        help="Prefix for generated output files",
        type=str,
        default="./velocity_"
    )
    return general_parser


def main(argsl=None):
    if argsl is None:
        argsl = sys.argv[1:]
    args, _ = arg_parser().parse_known_args(argsl)
    args = normalize_args(args, ["threads", "genes", "dpi", "top", "ndim", "nneighbors"])

    print("Setting scvelo parameters")
    configure_scvelo(args)

    print("Loading data")
    velocity_data = get_velocity_data(args)
    
    print("Estimating velocity")
    estimate_velocity(velocity_data, args)
    export_velocity_plot(velocity_data, args)
    export_text_data(velocity_data, args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))