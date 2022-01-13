#!/usr/bin/env python3
import os
import sys
import csv
import argparse
import yaml


def frac_to_perc(s):
    return round(float(s)*100.0, 1)


JOINT = {
    "Feature linkages detected": {
        "alias": "Feature linkages detected",
        "function": int
    },
    "Linked genes": {
        "alias": "Linked genes",
        "function": int
    },
    "Linked peaks": {
        "alias": "Linked peaks",
        "function": int
    },
    "Estimated number of cells": {
        "alias": "Estimated number of cells",
        "function": int
    },
    "order": [
        "Feature linkages detected",
        "Linked genes",
        "Linked peaks",
        "Estimated number of cells"
    ]
}

GEX_SEQUENCING = {
    "GEX Sequenced read pairs": {
        "alias": "Sequenced read pairs",
        "function": int
    },
    "GEX Valid barcodes": {
        "alias": "Valid barcodes, %",
        "function": frac_to_perc
    },
    "GEX Valid UMIs": {
        "alias": "Valid UMIs, %",
        "function": frac_to_perc
    },
    "GEX Q30 bases in barcode": {
        "alias": "Q30 bases in barcode, %",
        "function": frac_to_perc
    },
    "GEX Q30 bases in UMI": {
        "alias": "Q30 bases in UMI, %",
        "function": frac_to_perc
    },
    "GEX Q30 bases in read 2": {
        "alias": "Q30 bases in read 2, %",
        "function": frac_to_perc
    },
    "GEX Reads with TSO": {
        "alias": "Reads with TSO, %",
        "function": frac_to_perc
    },
    "GEX Percent duplicates": {
        "alias": "Percent duplicates, %",
        "function": frac_to_perc
    },
    "order": [
        "Sequenced read pairs",
        "Valid barcodes, %",
        "Valid UMIs, %",
        "Q30 bases in barcode, %",
        "Q30 bases in UMI, %",
        "Q30 bases in read 2, %",
        "Reads with TSO, %",
        "Percent duplicates, %"
    ]
}

GEX_CELLS = {
    "GEX Mean raw reads per cell": {
        "alias": "Mean raw reads per cell",
        "function": int
    },
    "GEX Fraction of transcriptomic reads in cells": {
        "alias": "Fraction of transcriptomic reads in cells, %",
        "function": frac_to_perc
    },
    "GEX Median UMI counts per cell": {
        "alias": "Median UMI counts per cell",
        "function": int
    },
    "GEX Median genes per cell": {
        "alias": "Median genes per cell",
        "function": int
    },
    "GEX Total genes detected": {
        "alias": "Total genes detected",
        "function": int
    },
    "order": [
        "Mean raw reads per cell",
        "Fraction of transcriptomic reads in cells, %",
        "Median UMI counts per cell",
        "Median genes per cell",
        "Total genes detected"
    ]
}

GEX_MAPPING = {
    "GEX Reads mapped to genome": {
        "alias": "Reads mapped to genome, %",
        "function": frac_to_perc
    },
    "GEX Reads mapped confidently to genome": {
        "alias": "Reads mapped confidently to genome, %",
        "function": frac_to_perc
    },
    "GEX Reads mapped confidently to intergenic regions": {
        "alias": "Reads mapped confidently to intergenic regions, %",
        "function": frac_to_perc
    },
    "GEX Reads mapped confidently to intronic regions": {
        "alias": "Reads mapped confidently to intronic regions, %",
        "function": frac_to_perc
    },
    "GEX Reads mapped confidently to exonic regions": {
        "alias": "Reads mapped confidently to exonic regions, %",
        "function": frac_to_perc
    },
    "GEX Reads mapped confidently to transcriptome": {
        "alias": "Reads mapped confidently to transcriptome, %",
        "function": frac_to_perc
    },
    "GEX Reads mapped antisense to gene": {
        "alias": "Reads mapped antisense to gene, %",
        "function": frac_to_perc
    },
    "order": [
        "Reads mapped to genome, %",
        "Reads mapped confidently to genome, %",
        "Reads mapped confidently to intergenic regions, %",
        "Reads mapped confidently to intronic regions, %",
        "Reads mapped confidently to exonic regions, %",
        "Reads mapped confidently to transcriptome, %",
        "Reads mapped antisense to gene, %"
    ]
}

ATAC_SEQUENCING = {
    "ATAC Sequenced read pairs": {
        "alias": "Sequenced read pairs",
        "function": int
    },
    "ATAC Valid barcodes": {
        "alias": "Valid barcodes, %",
        "function": frac_to_perc
    },
    "ATAC Q30 bases in barcode": {
        "alias": "Q30 bases in barcode, %",
        "function": frac_to_perc
    },
    "ATAC Q30 bases in read 1": {
        "alias": "Q30 bases in read 1, %",
        "function": frac_to_perc
    },
    "ATAC Q30 bases in read 2": {
        "alias": "Q30 bases in read 2, %",
        "function": frac_to_perc
    },
    "ATAC Percent duplicates": {
        "alias": "Percent duplicates, %",
        "function": frac_to_perc
    },
    "order": [
        "Sequenced read pairs",
        "Valid barcodes, %",
        "Q30 bases in barcode, %",
        "Q30 bases in read 1, %",
        "Q30 bases in read 2, %",
        "Percent duplicates, %"
    ]
}

ATAC_CELLS = {
    "ATAC Mean raw read pairs per cell": {
        "alias": "Mean raw read pairs per cell",
        "function": int
    },
    "ATAC Fraction of high-quality fragments in cells": {
        "alias": "Fraction of high-quality fragments in cells, %",
        "function": frac_to_perc
    },
    "ATAC Fraction of transposition events in peaks in cells": {
        "alias": "Fraction of transposition events in peaks in cells, %",
        "function": frac_to_perc
    },
    "ATAC Median high-quality fragments per cell": {
        "alias": "Median high-quality fragments per cell",
        "function": int
    },
    "order": [
        "Mean raw read pairs per cell",
        "Fraction of high-quality fragments in cells, %",
        "Fraction of transposition events in peaks in cells, %",
        "Median high-quality fragments per cell"
    ]
}

ATAC_TARGETING = {
    "ATAC Number of peaks": {
        "alias": "Number of peaks",
        "function": int
    },
    "ATAC Fraction of genome in peaks": {
        "alias": "Fraction of genome in peaks, %",
        "function": frac_to_perc
    },
    "ATAC TSS enrichment score": {
        "alias": "TSS enrichment score",
        "function": float
    },
    "ATAC Fraction of high-quality fragments overlapping TSS": {
        "alias": "Fraction of high-quality fragments overlapping TSS, %",
        "function": frac_to_perc
    },
    "ATAC Fraction of high-quality fragments overlapping peaks": {
        "alias": "Fraction of high-quality fragments overlapping peaks, %",
        "function": frac_to_perc
    },
    "order": [
        "Number of peaks",
        "Fraction of genome in peaks, %",
        "TSS enrichment score",
        "Fraction of high-quality fragments overlapping TSS, %",
        "Fraction of high-quality fragments overlapping peaks, %"
    ]
}

ATAC_MAPPING = {
    "ATAC Confidently mapped read pairs": {
        "alias": "Confidently mapped read pairs, %",
        "function": frac_to_perc
    },
    "ATAC Unmapped read pairs": {
        "alias": "Unmapped read pairs, %",
        "function": frac_to_perc
    },
    "ATAC Non-nuclear read pairs": {
        "alias": "Non-nuclear read pairs, %",
        "function": frac_to_perc
    },
    "order": [
        "Confidently mapped read pairs, %",
        "Unmapped read pairs, %",
        "Non-nuclear read pairs, %"
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
    headers = ["General", "GEX Sequencing", "GEX Cells", "GEX Mapping", "ATAC Sequencing", "ATAC Cells", "ATAC Targeting", "ATAC Mapping"]
    ref_dicts = [JOINT, GEX_SEQUENCING, GEX_CELLS, GEX_MAPPING, ATAC_SEQUENCING, ATAC_CELLS, ATAC_TARGETING, ATAC_MAPPING]
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