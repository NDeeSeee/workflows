import os
import sys
import argparse
import anndata
import scvelo
import pandas
import numpy
import matplotlib


def estimate_velocity(args):
    velocity_data = anndata.read_h5ad(args.h5ad)
    print(velocity_data)
    cells_data = pandas.read_csv(args.cells, header=None, sep="\t")
    umap_data = pandas.read_csv(args.umap, header=None, sep="\t")
    clusters_data = pandas.read_csv(args.clusters, header=None, sep="\t")
    filtered_velocity_data = velocity_data[cells_data.iloc[:, 0]].copy()
    print(filtered_velocity_data)
    filtered_velocity_data.obsm["X_umap"] = umap_data.values
    filtered_velocity_data.obs["cluster"] = clusters_data.values
    scvelo.pl.proportions(
        filtered_velocity_data,
        show=False,
        save=".png"
    )
    print(filtered_velocity_data)
    scvelo.pp.filter_and_normalize(filtered_velocity_data)
    scvelo.pp.moments(filtered_velocity_data)
    scvelo.tl.velocity(filtered_velocity_data, mode="stochastic")
    scvelo.tl.velocity_graph(filtered_velocity_data, n_jobs=args.threads)
    scvelo.pl.velocity_embedding_stream(
        filtered_velocity_data,
        basis="umap",
        color="cluster",
        show=False,
        save="stream.png"
    )
    scvelo.pl.velocity_embedding_grid(
        filtered_velocity_data,
        basis="umap",
        color="cluster",
        show=False,
        save="grid.png"
    )
    if len(args.genes) > 0:
        scvelo.pl.velocity(
            filtered_velocity_data,
            args.genes,
            show=False,
            save="genes.png"
        )
    scvelo.tl.velocity_confidence(filtered_velocity_data)
    scvelo.pl.scatter(
        filtered_velocity_data,
        c=["velocity_length", "velocity_confidence"],
        cmap="coolwarm",
        perc=[5, 95],
        show=False,
        save="confidence.png"
    )
    scvelo.tl.velocity_pseudotime(filtered_velocity_data)
    scvelo.pl.scatter(
        filtered_velocity_data,
        color="velocity_pseudotime",
        cmap="gnuplot",
        show=False,
        save="pseudotime.png"
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
    general_parser.add_argument("--h5ad",     help="Path to the h5ad file",                    type=str, required=True)
    general_parser.add_argument("--cells",    help="Path to the cells data TSV file",          type=str, required=True)
    general_parser.add_argument("--umap",     help="Path to the UMAP data TSV file",           type=str, required=True)
    general_parser.add_argument("--genes",    help="List of genes of interest",                type=str, nargs="*")
    general_parser.add_argument("--clusters", help="Path to the clusters data TSV file",       type=str, required=True)
    general_parser.add_argument("--threads",  help="Number of threads to run in parallel",     type=int, default=1)
    return general_parser


def main(argsl=None):
    if argsl is None:
        argsl = sys.argv[1:]
    args, _ = arg_parser().parse_known_args(argsl)
    args = normalize_args(args, ["threads", "genes"])
    scvelo.set_figure_params("scvelo", dpi=300, transparent=False)
    estimate_velocity(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))