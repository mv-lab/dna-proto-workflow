#!/usr/bin/env python3
from sys import stderr, stdout, stdin, argv
from os.path import basename

def extractSN(filename):
    fields = {}
    with open(filename) as fh:
        for line in fh:
            if not line.startswith("SN"):
                continue
            row = line.rstrip('\n').split('\t')
            field = row[1].rstrip(":").replace("(", ""). \
                           replace(")", "").replace(" ", "_")
            val = row[2]
            fields[field] = val
    return fields

def main():
    fields = None
    for file in argv[1:]:
        f = extractSN(file)
        if fields is None:
            fields = ["sample"] + list(f.keys())
            print(*fields, sep="\t")

        row = [basename(file), ]
        for field in fields[1:]: # skip sample
            row.append(f.get(field, "NA"))
        print(*row, sep="\t")


if __name__ == "__main__":
    main()
