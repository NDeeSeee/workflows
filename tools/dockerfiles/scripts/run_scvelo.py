import os
import sys
import argparse
import anndata
import scvelo
import pandas
import numpy
import matplotlib
import dataframe_image


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


def get_velocity_data(h5ad, cells, umap, pca, clusters, barcode_suffix=None):
    raw_velocity_data = anndata.read_h5ad(h5ad)
    cells_data = pandas.read_csv(cells, header=None, sep="\t")
    filtered_velocity_data = raw_velocity_data[cells_data.iloc[:, 0]].copy()
    clusters_data = pandas.read_csv(clusters, header=None, sep="\t")
    filtered_velocity_data.obs["clusters"] = clusters_data.values
    # https://github.com/theislab/scvelo/issues/37#issuecomment-499101214
    # https://github.com/theislab/scvelo/discussions/722
    pca_data = pandas.read_csv(pca, header=None, sep="\t")
    umap_data = pandas.read_csv(umap, header=None, sep="\t")
    filtered_velocity_data.obsm["X_pca"] = pca_data.values
    filtered_velocity_data.obsm["X_umap"] = umap_data.values
    if barcode_suffix != None:
        filtered_velocity_data.obs.index = filtered_velocity_data.obs.index + "-" + barcode_suffix
    return filtered_velocity_data


def estimate_velocity(velocity_data, args):
    scvelo.pp.filter_and_normalize(velocity_data)
    scvelo.pp.moments(velocity_data, n_pcs=None, n_neighbors=args.nneighbors)
    scvelo.tl.recover_dynamics(velocity_data, n_jobs=args.threads)
    scvelo.tl.velocity(velocity_data, mode="dynamical")
    scvelo.tl.velocity_graph(velocity_data, n_jobs=args.threads)
    scvelo.tl.velocity_confidence(velocity_data)
    scvelo.tl.rank_dynamical_genes(velocity_data, groupby="clusters", n_genes=args.top)
    scvelo.tl.latent_time(velocity_data)
    velocity_data.uns["neighbors"]["distances"] = velocity_data.obsp["distances"]
    velocity_data.uns["neighbors"]["connectivities"] = velocity_data.obsp["connectivities"]
    scvelo.tl.paga(
        velocity_data,
        groups="clusters",
        use_time_prior="latent_time"
    )


def export_h5ad_data(velocity_data, args):
    velocity_data.write(
        filename=args.output + "compressed.h5ad",
        compression="gzip"
    )


def export_text_data(velocity_data, args):
    driver_genes = scvelo.get_df(velocity_data, "rank_dynamical_genes/names")
    driver_genes.to_csv(args.output + "putative_driver_genes.tsv", sep="\t")


def export_velocity_plot(velocity_data, args):
    try:
        scvelo.pl.proportions(
           velocity_data,
            groupby="clusters",
            show=False,              # need to specify explicitely for scvelo.pl.proportions
            dpi=args.dpi,            # need to specify explicitely for scvelo.pl.proportions
            save="chart.png"         # "proportions" part is hardcoded so we add only extension
        )
    except Exception as err:
        print("Failed to export spliced/unspliced proportions charts\n", err)
    
    # These provide insights where cells differentiate at a slower/faster pace, and where the direction is un-/determined.
    try:
        scvelo.pl.scatter(
            velocity_data,
            basis="umap",
            color="velocity_length",
            color_map="coolwarm",
            colorbar=True,
            fontsize=10,
            legend_fontsize=8,
            title="Velocity length",
            add_text="defines speed of differentiation",
            add_text_pos=(0.05, 0),
            save="velocity_length.png"
        )
    except Exception as err:
        print("Failed to export velocity length plot\n", err)
    
    try:
        scvelo.pl.scatter(
            velocity_data,
            basis="umap",
            color="velocity_confidence",
            color_map="coolwarm",
            colorbar=True,
            fontsize=10,
            legend_fontsize=8,
            title="Velocity confidence",
            add_text="defines velocity vectors direction correlation",
            add_text_pos=(0.05, 0),
            save="velocity_confidence.png"
        )
    except Exception as err:
        print("Failed to export velocity confidence plot\n", err)

    try:
        velocity_speed_df = velocity_data.obs.groupby("clusters")["velocity_length", "velocity_confidence"].mean().T
        dataframe_image.export(
            velocity_speed_df.style.background_gradient(cmap="coolwarm", axis=1),
            args.output + "velocity_metrics.png"
        )
    except Exception as err:
        print("Failed to export velocity length and confidence metrics\n", err)

    # These yields fine-grained insights into the developmental processes
    try:
        scvelo.pl.velocity_embedding_stream(
            velocity_data,
            basis="umap",
            color="clusters",
            fontsize=10,
            legend_fontsize=8,
            title="Velocities streams",
            save="velocity_stream.png"
        )
    except Exception as err:
        print("Failed to export velocities stream plot\n", err)

    try:
        scvelo.pl.velocity_embedding_grid(
            velocity_data,
            basis="umap",
            color="clusters",
            fontsize=10,
            legend_fontsize=8,
            title="Velocities grid",
            save="velocity_grid.png"
        )
    except Exception as err:
        print("Failed to export velocities grid plot\n", err)

    try:
        scvelo.pl.scatter(
            velocity_data,
            basis="umap",
            color="latent_time",
            color_map="gnuplot",
            alpha=0.25,
            fontsize=10,
            legend_fontsize=8,
            title="Latent time",
            add_text="defines the real time experienced by cells as they differentiate\n(based on the transcriptional dynamics)",
            add_text_pos=(0.05, 0),
            save="latent_time.png"
        )
    except Exception as err:
        print("Failed to export latent time plot\n", err)

    top_genes = velocity_data.var["fit_likelihood"].sort_values(ascending=False).index[:args.top]
    try:
        scvelo.pl.heatmap(
            velocity_data,
            var_names=top_genes,
            sortby="latent_time",
            n_convolve=30,
            col_color="clusters",
            font_scale=0.2,
            colorbar=True,
            yticklabels=True,
            save="driver_genes.png"
        )
    except Exception as err:
        print("Failed to export driver genes heatmap\n", err)

    for gene in top_genes:
        try:
            scvelo.pl.velocity(
                velocity_data,
                var_names=gene,
                basis="umap",
                color="clusters",
                dpi=args.dpi,
                save=f"""phase_{gene}.png"""
            )
        except Exception as err:
            print("Failed to export driver gene phase portraits\n", err)
        try:
            scvelo.pl.scatter(
                velocity_data,
                basis="umap",
                color="clusters",
                x="latent_time",
                y=gene,
                frameon=False,
                fontsize=10,
                legend_fontsize=8,
                legend_loc="best",
                title=f"""{gene} expression dynamics""",
                dpi=args.dpi,
                save=f"""expression_{gene}.png"""
            )
        except Exception as err:
            print("Failed to export driver gene expression dynamics plot", err)

    try:
        scvelo.pl.paga(
            velocity_data,
            basis="umap",
            alpha=0.1,
            min_edge_width=2,
            node_size_scale=1.5,
            size=50,
            legend_loc="on data",
            fontsize=10,
            legend_fontsize=8,
            fontoutline=1,
            title="PAGA graph with velocity-directed edges",
            add_text="abstracted graph of partitions, in which edge weights represent\nconfidence in the presence of connections",
            add_text_pos=(0.05, 0),
            save="paga.png"
        )
    except Exception as err:
        print("Failed to export PAGA plot", err)

    # if user provided genes of interest
    for gene in args.genes or []:
        try:
            scvelo.pl.velocity(
                velocity_data,
                var_names=gene,
                basis="umap",
                color="clusters",
                dpi=args.dpi,
                save=f"""custom_phase_{gene}.png"""
            )
        except Exception as err:
            print("Failed to export driver gene phase portraits\n", err)
        try:
            scvelo.pl.scatter(
                velocity_data,
                basis="umap",
                color="clusters",
                x="latent_time",
                y=gene,
                frameon=False,
                fontsize=10,
                legend_fontsize=8,
                legend_loc="best",
                title=f"""{gene} expression dynamics""",
                dpi=args.dpi,
                save=f"""custom_expression_{gene}.png"""
            )
        except Exception as err:
            print("Failed to export driver gene expression dynamics plot", err)


def normalize_args(args, skip_list=[]):
    normalized_args = {}
    for key,value in args.__dict__.items():
        if key not in skip_list:
            if isinstance(value, list):
                normalized_args[key] = [
                    item if not item or os.path.isabs(item) else os.path.normpath(os.path.join(os.getcwd(), item))
                    for item in value
                ]
            else:
                normalized_args[key] = value if not value or os.path.isabs(value) else os.path.normpath(os.path.join(os.getcwd(), value))
        else:
            normalized_args[key]=value
    return argparse.Namespace (**normalized_args)


def arg_parser():
    general_parser = argparse.ArgumentParser()
    general_parser.add_argument(
        "--h5ad",
        help="Path to the h5ad file. \
              If multiple files provided, the order should correspond to \
              the order of files provided in --cells, --umap, --pca, --clusters",
        type=str,
        required=True,
        nargs="*"
    )
    general_parser.add_argument(
        "--cells",
        help="Path to the cells data TSV file. \
              If multiple files provided, the order should correspond to \
              the order of files provided in --h5ad, --umap, --pca, --clusters",
        type=str,
        required=True,
        nargs="*"
    )
    general_parser.add_argument(
        "--umap",
        help="Path to the custom UMAP data TSV file. \
              If multiple files provided, the order should correspond to \
              the order of files provided in --cells, --h5ad, --pca, --clusters",
        type=str,
        required=True,
        nargs="*"
    )
    general_parser.add_argument(
        "--pca",
        help="Path to the custom PCA data TSV file.\
              If multiple files provided, the order should correspond to \
              the order of files provided in --cells, --umap, --h5ad, --clusters",
        type=str,
        required=True,
        nargs="*"
    )
    general_parser.add_argument(
        "--clusters",
        help="Path to the clusters data TSV file. \
              If multiple files provided, the order should correspond to \
              the order of files provided in --cells, --umap, --pca, --h5ad",
        type=str,
        required=True,
        nargs="*"
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
    args = normalize_args(args, ["threads", "genes", "dpi", "top", "nneighbors"])

    print("Setting scvelo parameters")
    configure_scvelo(args)

    concat_velocity_data = []                                                                                     # can't use None, because anndata can't be compred to None. We'll use len() == 0 instead
    for barcode_suffix, (h5ad, cells, umap, pca, clusters) in enumerate(zip(args.h5ad, args.cells, args.umap, args.pca, args.clusters)):
        barcode_suffix = str(barcode_suffix + 1) if len(args.h5ad) > 1 else None                                  # safety measure in case we don't need to concatenate input h5ad files
        print("Loading dataset")
        print(f"""  velocity - {h5ad}\n  cells    - {cells}\n  UMAP     - {umap}\n  PCA      - {pca}\n  clusters - {clusters}""")
        velocity_data = get_velocity_data(h5ad, cells, umap, pca, clusters, barcode_suffix)
        if len(concat_velocity_data) == 0:
            concat_velocity_data = velocity_data
        else:
            concat_velocity_data = concat_velocity_data.concatenate(
                velocity_data,
                batch_key="dataset",
                index_unique=None
            )

    print("Estimating velocity")
    estimate_velocity(concat_velocity_data, args)
    export_h5ad_data(concat_velocity_data, args)
    export_text_data(concat_velocity_data, args)
    export_velocity_plot(concat_velocity_data, args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))