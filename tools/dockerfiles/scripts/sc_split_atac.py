#!/usr/bin/env python
import os
import gc
import sys
import time
import logging
import argparse
import numpy as np
import pandas as pd
import multiprocessing
from functools import partial
from collections import namedtuple
from multiprocessing.shared_memory import SharedMemory
from multiprocessing.managers import SharedMemoryManager


Job = namedtuple("Job", "cluster sh_mem_block_name np_data_shape np_data_dtype")


def get_fragments_and_clusters(fragments_location, clusters_location):
    logging.info("Loading fragments data")
    fragments_data = pd.read_csv(
        fragments_location,
        sep="\t", 
        names=["chr", "start", "end", "barcode", "reads"],
        compression="gzip"
    )
    logging.debug(f"\n{fragments_data.head()}")
    logging.info("Loading clusters data")
    clusters_data = pd.read_csv(
        clusters_location,
        sep="\t",
        names=["barcode", "cluster"],
    )
    logging.debug(f"\n{clusters_data.head()}")
    logging.info("Aggregating fragments and clusters")
    aggregated_data = pd.merge(                                            # inner join
        left=fragments_data,
        right=clusters_data,
        left_on="barcode",
        right_on="barcode"
    )
    logging.debug(f"\n{aggregated_data.head()}")
    longest_chr = aggregated_data["chr"].map(lambda x: len(str(x))).max()
    longest_barcode = aggregated_data["barcode"].map(lambda x: len(str(x))).max()
    longest_cluster = aggregated_data["cluster"].map(lambda x: len(str(x))).max()
    logging.debug("Converting aggregated data to NumPy RecArray")
    np_data = aggregated_data.to_records(
        index=False,
        column_dtypes={                                                    # need to explicitly set types, otherwise causes segmenation fault
            "chr": f"U{longest_chr}",
            "barcode": f"U{longest_barcode}",
            "cluster": f"U{longest_cluster}"
        }
    )
    logging.debug(f"  size:  {np_data.nbytes/1e6}MB")
    logging.debug(f"  shape: {np_data.shape}")
    logging.debug(f"  dtype: {np_data.dtype}")
    del aggregated_data
    del fragments_data
    del clusters_data
    gc.collect()
    return np_data


def process_cluster(args, job, np_data=None):
    setup_logger(multiprocessing.get_logger(), args.loglevel)
    multiprocessing.current_process().name = f"job_{job.cluster}"
    output_location = f"{args.output}_cluster_{job.cluster}.bed"
    logging.info(f"Extracting transcripts from {job.cluster} cluster to {output_location}")
    try:
        if np_data is None:
            shared_memory_block = SharedMemory(job.sh_mem_block_name)
            logging.debug(f"Accessing shared memory block {shared_memory_block.name}")
            np_data = np.recarray(
                shape=job.np_data_shape,
                dtype=job.np_data_dtype,
                buf=shared_memory_block.buf
            )
        selected_trascripts = pd.DataFrame(np_data[np_data["cluster"]==job.cluster])
        selected_trascripts[["chr", "start", "end", "barcode", "reads"]].sort_values(["chr", "start", "end"]).to_csv(
            output_location,
            sep="\t",
            header=False,
            index=False
        )
        del selected_trascripts
        gc.collect()
    except Exception as err:
        logging.warning(f"Failed to export {output_location} because of {err}")


def setup_logger(logger, log_level, log_format=None):
    log_format = "%(processName)12s (%(asctime)s): %(message)s" if log_format is None else log_format
    for log_handler in logger.handlers:
        logger.removeHandler(log_handler)
    for log_filter in logger.filters:
        logger.removeFilter(log_filter)
    logging.basicConfig(level=log_level, format=log_format)


def get_jobs(args, sh_mem_block_name, np_data_shape, np_data_dtype):
    logging.info("Preparing jobs")
    clusters_data = pd.read_csv(
        args.clusters,
        sep="\t",
        names=["barcode", "cluster"]
    )
    return [
        Job(
            cluster=str(cluster),
            sh_mem_block_name=sh_mem_block_name,
            np_data_shape=np_data_shape,
            np_data_dtype=np_data_dtype
        )
        for cluster in sorted(clusters_data["cluster"].unique())
    ]


class ArgsParser():

    def __init__(self, args):
        self.args, _ = self.get_parser().parse_known_args(args)
        self.normalize_args(["cpus", "loglevel"])
        self.assert_args()
        self.set_args_as_attributes()

    def set_args_as_attributes(self):
        for arg, value in vars(self.args).items():
            setattr(self, arg, value)

    def get_parser(self):
        general_parser = argparse.ArgumentParser()
        general_parser.add_argument(
            "--fragments",
            help="Path to GZIP compressed TSV file with ATAC fragments (from Cell Ranger ARC)",
            type=str, required=True
        )
        general_parser.add_argument(
            "--clusters",
            help="Path to headerless TSV file with barcodes (first column) and clusters (second column)",
            type=str, required=True
        )
        general_parser.add_argument(
            "--cpus",
            help=" ".join(
                [
                    "Number of processes to run in parallel. When set to more than 1,",
                    "uses shared memory for interprocess communication.",
                    "Default: 1 (do not use shared memory)"
                ]
            ),
            type=int, default=1
        )
        general_parser.add_argument(
            "--loglevel",
            help="Logging level. Default: info",
            type=str, default="info",
            choices=["fatal", "error", "warning", "info", "debug"]
        )
        general_parser.add_argument(
            "--output",
            help="Output file prefix",
            type=str, default="split"
        )
        return general_parser

    def normalize_args(self, skip=None):
        skip = [] if skip is None else skip
        normalized_args = {}
        for key,value in self.args.__dict__.items():
            if key not in skip:
                normalized_args[key] = value if not value or os.path.isabs(value) else os.path.normpath(os.path.join(os.getcwd(), value))
            else:
                normalized_args[key]=value
        self.args = argparse.Namespace (**normalized_args)

    def assert_args(self):
        self.args.loglevel = {
            "fatal": logging.FATAL,
            "error": logging.ERROR,
            "warning": logging.WARNING,
            "info": logging.INFO,
            "debug": logging.DEBUG
        }[self.args.loglevel]


def main(args=None):
    args = ArgsParser(sys.argv[1:] if args is None else args)
    setup_logger(logging.root, args.loglevel)
    start = time.time()
    np_data = get_fragments_and_clusters(args.fragments, args.clusters)

    if args.cpus == 1:
        jobs = get_jobs(args, None, None, None)                                                   # no need to set parameters if SharedMemory not used
        logging.info(f"Running {len(jobs)} job(s) using {args.cpus} CPU (no shared memory)")
        for job in jobs:
            process_cluster(args, job, np_data)
    else:
        with SharedMemoryManager() as sm_manager:
            shared_memory_block = sm_manager.SharedMemory(size=np_data.nbytes)                    # create a new shared memory block
            logging.debug(f"Setting up shared memory block {shared_memory_block.name}")
            shared_np_data = np.recarray(                                                         # create a NumPy array backed by shared memory
                shape=np_data.shape,
                dtype=np_data.dtype,
                buf=shared_memory_block.buf
            )
            logging.debug(f"Copying data to the shared memory block {shared_memory_block.name}")
            np.copyto(shared_np_data, np_data)                                                    # copy the original data into shared memory block
            jobs = get_jobs(args, shared_memory_block.name, np_data.shape, np_data.dtype)
            del np_data
            gc.collect()
            logging.info(f"Running {len(jobs)} job(s) using {args.cpus} CPUs (with shared memory)")
            with multiprocessing.Pool(args.cpus) as pool:
                pool.map(partial(process_cluster, args), jobs)
    logging.info(f"Elapsed time: {round(time.time() - start)} sec")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))