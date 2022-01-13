#!/usr/bin/env python3
import os
import sys
import csv
import argparse
import yaml


def cut_percent(s):
    return round(float(s.strip().replace("%", "").split()[0]), 1)

def clean_int(s):
    return int(s.replace(",", ""))


GEX_SEQUENCING = {
    "Number of Reads": {
        "alias": "Number of Reads",
        "function": clean_int
    },
    "Valid Barcodes": {
        "alias": "Valid Barcodes, %",
        "function": cut_percent
    },
    "Sequencing Saturation": {
        "alias": "Sequencing Saturation, %",
        "function": cut_percent
    },
    "Q30 Bases in Barcode": {
        "alias": "Q30 Bases in Barcode, %",
        "function": cut_percent
    },
    "Q30 Bases in RNA Read": {
        "alias": "Q30 Bases in RNA Read, %",
        "function": cut_percent
    },
    "Q30 Bases in UMI": {
        "alias": "Q30 Bases in UMI, %",
        "function": cut_percent
    },
    "order": [
        "Number of Reads",
        "Valid Barcodes, %",
        "Sequencing Saturation, %",
        "Q30 Bases in Barcode, %",
        "Q30 Bases in RNA Read, %",
        "Q30 Bases in UMI, %"
    ]
}

GEX_CELLS = {
    "Estimated Number of Cells": {
        "alias": "Estimated Number of Cells",
        "function": clean_int
    },
    "Fraction Reads in Cells": {
        "alias": "Fraction Reads in Cells, %",
        "function": cut_percent
    },
    "Mean Reads per Cell": {
        "alias": "Mean Reads per Cell",
        "function": clean_int
    },
    "Median Genes per Cell": {
        "alias": "Median Genes per Cell",
        "function": clean_int
    },
    "Total Genes Detected": {
        "alias": "Total Genes Detected",
        "function": clean_int
    },
    "Median UMI Counts per Cell": {
        "alias": "Median UMI Counts per Cell",
        "function": clean_int
    },
    "order": [
        "Estimated Number of Cells",
        "Fraction Reads in Cells, %",
        "Mean Reads per Cell",
        "Median Genes per Cell",
        "Total Genes Detected",
        "Median UMI Counts per Cell"
    ]
}

GEX_MAPPING = {
    "Reads Mapped to Genome": {
        "alias": "Reads Mapped to Genome, %",
        "function": cut_percent
    },
    "Reads Mapped Confidently to Genome": {
        "alias": "Reads Mapped Confidently to Genome, %",
        "function": cut_percent
    },
    "Reads Mapped Confidently to Intergenic Regions": {
        "alias": "Reads Mapped Confidently to Intergenic Regions, %",
        "function": cut_percent
    },
    "Reads Mapped Confidently to Intronic Regions": {
        "alias": "Reads Mapped Confidently to Intronic Regions, %",
        "function": cut_percent
    },
    "Reads Mapped Confidently to Exonic Regions": {
        "alias": "Reads Mapped Confidently to Exonic Regions, %",
        "function": cut_percent
    },
    "Reads Mapped Confidently to Transcriptome": {
        "alias": "Reads Mapped Confidently to Transcriptome, %",
        "function": cut_percent
    },
    "Reads Mapped Antisense to Gene": {
        "alias": "Reads Mapped Antisense to Gene, %",
        "function": cut_percent
    },
    "order": [
        "Reads Mapped to Genome, %",
        "Reads Mapped Confidently to Genome, %",
        "Reads Mapped Confidently to Intergenic Regions, %",
        "Reads Mapped Confidently to Intronic Regions, %",
        "Reads Mapped Confidently to Exonic Regions, %",
        "Reads Mapped Confidently to Transcriptome, %",
        "Reads Mapped Antisense to Gene, %"
    ]
}


def open_file(location):
    lines = []
    with open(location, "r") as input_stream:
        for line in input_stream:
            if line.strip():
                lines.append(line.strip())
    return lines


def get_key(data, long_k):
    for short_k, v in data.items():
        if short_k in long_k:
            return v["alias"], v["function"]
    raise Exception


def process_cell_ranger_report(location, collected_results, header, ref_dict):
    if not collected_results.get(header, None):
        collected_results[header] = {}
    keys, values = None, None
    with open(location, "r") as input_stream:
        for i, row in enumerate(csv.reader(input_stream)):
            if i==0:
                keys = row
            else:
                values = row
    for key, value in zip(keys, values):
        try:
            ref_key, ref_function = get_key(ref_dict, key)
            if not collected_results[header].get(ref_key, None):
                collected_results[header][ref_key] = ref_function(value)
        except Exception:
            pass            
    if ref_dict.get("order", None):
        collected_results[header] = {k: collected_results[header][k] for k in ref_dict["order"] if k in collected_results[header]}


def collect_stats(args):
    collected_results = {}
    headers = ["GEX Sequencing", "GEX Cells", "GEX Mapping"]
    ref_dicts = [GEX_SEQUENCING, GEX_CELLS, GEX_MAPPING]
    for header, ref_dict in zip(headers, ref_dicts):
        process_cell_ranger_report(
            location=args.metrics,
            collected_results=collected_results,
            header=header,
            ref_dict=ref_dict
        )
    return (collected_results)


def export_results_yaml(collected_data, location):
    with open(location+".yaml", "w") as output_stream:
        output_stream.write(yaml.dump(collected_data, width=1000, sort_keys=False))


def export_results_markdown(collected_data, location):
    with open(location+".md", "w") as output_stream:
        for line in yaml.dump(collected_data, width=1000, sort_keys=False).split("\n"):
            if not line.strip():
                continue
            if line.startswith("  - "):
                output_stream.write(line+"\n")
            elif line.startswith("    "):
                output_stream.write("<br>"+line+"\n")
            elif line.startswith("  "):
                output_stream.write("- "+line+"\n")
            else:
                output_stream.write("### "+line+"\n")


def export_results_table(collected_data, location):
    with open(location+".tsv", "w") as output_stream:
        header = []
        data = []
        for group_key, group_value in collected_data.items():
            header.append(group_key)
            data.append("")
            for key, value in group_value.items():
                header.append(key)
                data.append(value)
        output_stream.write("\t".join(header)+"\n")
        output_stream.write("\t".join(str(l) for l in data))


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
    general_parser.add_argument("--metrics", help="Path to the Cell Ranger (ARC) Count summary.csv file", required=True)
    general_parser.add_argument("--output",  help="Output filename prefix", default="collected_stats")
    return general_parser


def main(argsl=None):
    if argsl is None:
        argsl = sys.argv[1:]
    args,_ = arg_parser().parse_known_args(argsl)
    args = normalize_args(args, skip_list=["paired"])
    collected_data = collect_stats(args)
    export_results_yaml(collected_data, args.output)
    export_results_table(collected_data, args.output)
    export_results_markdown(collected_data, args.output)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))