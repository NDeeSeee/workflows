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


def get_filtered_velocity_data(velocity_data, cells_data, umap_data, clusters_data):
    filtered_velocity_data = velocity_data[cells_data.iloc[:, 0]].copy()
    filtered_velocity_data.obsm["X_umap"] = umap_data.values
    filtered_velocity_data.obs["cluster"] = clusters_data.values
    return filtered_velocity_data


def estimate_velocity(velocity_data, args):
    scvelo.pp.filter_and_normalize(velocity_data)
    scvelo.pp.moments(velocity_data)
    if args.mode == "dynamical":
        scvelo.tl.recover_dynamics(velocity_data, n_jobs=args.threads)
    scvelo.tl.velocity(velocity_data, mode=args.mode)
    scvelo.tl.velocity_graph(velocity_data, n_jobs=args.threads)
    scvelo.tl.velocity_confidence(velocity_data)
    scvelo.tl.velocity_pseudotime(velocity_data)


def export_velocity_plot(velocity_data, args):
    scvelo.pl.velocity_embedding_stream(
        velocity_data,
        basis="umap",
        color="cluster",
        save="stream.png"
    )
    scvelo.pl.velocity_embedding_grid(
        velocity_data,
        basis="umap",
        color="cluster",
        save="grid.png"
    )
    if args.genes and len(args.genes) > 0:
        scvelo.pl.velocity(
            velocity_data,
            args.genes,
            save="genes.png"
        )
    scvelo.pl.scatter(
        velocity_data,
        c=["velocity_length", "velocity_confidence"],
        cmap="coolwarm",
        perc=[5, 95],
        save="confidence.png"
    )
    scvelo.pl.scatter(
        velocity_data,
        color="velocity_pseudotime",
        cmap="gnuplot",
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
        help="Path to the UMAP data TSV file",
        type=str,
        required=True
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
        "--mode",
        help="Velocity calculation mode",
        type=str,
        choices=["stochastic", "deterministic", "dynamical"],
        default="stochastic"
    )
    general_parser.add_argument(
        "--threads",
        help="Number of threads to run in parallel",
        type=int,
        default=1
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
    args = normalize_args(args, ["threads", "genes", "dpi", "mode"])
    configure_scvelo(args)

    print(f"""Loading data:\n velocity - {args.h5ad}\n barcodes - {args.cells}\n UMAP     - {args.umap}\n clusters - {args.clusters}""")
    raw_velocity_data = anndata.read_h5ad(args.h5ad)
    cells_data = pandas.read_csv(args.cells, header=None, sep="\t")
    umap_data = pandas.read_csv(args.umap, header=None, sep="\t")
    clusters_data = pandas.read_csv(args.clusters, header=None, sep="\t")

    filtered_velocity_data = get_filtered_velocity_data(raw_velocity_data, cells_data, umap_data, clusters_data)
    estimate_velocity(filtered_velocity_data, args)
    export_velocity_plot(filtered_velocity_data, args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))