#!/usr/bin/env python3
from sys import stderr, stdout, stdin, argv
from os.path import basename, splitext
import re
import json
import csv
import argparse

try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda x: x


def doSN(line, data):
    if "SN" not in data:
        data["SN"] = {}
    row = line.rstrip('\n').split('\t')
    field = re.sub(r'\(|\)|:|%', '', row[1])
    field = re.sub(r'\s+$|^\s+', '', field)
    field = re.sub(r'\s+', '_', field)
    val = row[2]
    data["SN"][field] = val


# TODO: not handling this right now
# def doTableByCycle(line, data, key):
#     cols = {
#         "FFQ": ["cycle_number", "quality"],
#         "LFQ": ["cycle_number", "quality"],
#     }


TABLE_COLS = {
    "IS": ["insert_size", "pairs", "pairs_in", "pairs_out", "pairs_other"],
    "ID": ["indel_size", "n_insertion", "n_deletion"],
    "GCF": ["gc_content", "n_reads"],
    "GCL": ["gc_content", "n_reads"],
    "GCC": ["cycle_number", "a_pct", "c_pct", "g_pct", "t_pct", "n_pct", None],
    "RL": ["read_length", "n_reads"],
    "COV": ["coverage_bin", None, "genome_bp"],
}

def doTable(line, data):
    key = line.strip().split('\t')[0]
    if key not in data:
        data[key] = []
    vals = line.strip().split("\t")
    valdict = {}
    for i, col in enumerate(TABLE_COLS[key]):
        if col is None:
            continue
        valdict[col] = vals[i+1]
    data[key].append(valdict)


def dump(samples, key, outfile):
    with open(outfile, "w") as fh:
        csvfh = None
        for sample, alldata in samples.items():
            if key not in alldata:
                continue
            data = alldata[key]
            if isinstance(data, dict):
                # for some where rows are split over multiple lines, e.g. SN
                data = [data, ]
            if csvfh is None:
                header = ["sample", ] + list(data[0].keys())
                csvfh = csv.DictWriter(fh, header)
                csvfh.writeheader()
            for row in data:
                row["sample"] = sample
                csvfh.writerow(row)



def main():
    p = argparse.ArgumentParser(prog="tidybamstat")
    p.add_argument("-o", "--outprefix", type=str, required=True,
                   help="Output prefix")
    p.add_argument("input", type=str, nargs='+')
    args = p.parse_args()

    samples = {}
    print("Load alignment stats:", file=stderr)
    for file in tqdm(args.input):
        data = dict()
        with open(file) as fh:
            for line in fh:
                dtype = line.strip().split("\t")[0]
                if line.startswith("#"):
                    continue
                elif line.startswith("SN"):
                    doSN(line, data)
                elif dtype in TABLE_COLS:
                    doTable(line, data)
        sample = splitext(basename(file))[0]
        samples[sample] = data

    allkeys = ["SN"] + list(TABLE_COLS.keys())
    print("Dumping data by data type: ", file=stderr, end='', flush=True)
    for key in allkeys:
        print(f"{key}, ", file=stderr, end='',  flush=True)
        outfile = f"{args.outprefix}_{key}.csv"
        dump(samples, key, outfile)
    print(" Done!", file=stderr)

if __name__  == "__main__":
    main()
