#!/usr/bin/env python3
import numpy as np
import argparse as ap

def main(src, dst):
    mat = np.load(src)
    if mat.shape[1] > mat.shape[0]:
        # Ensure we have a "long" table or we'll kill R
        mat = mat.T
    np.savetxt(dst, mat, delimiter="\t")

if __name__ == "__main__":
    a = ap.ArgumentParser("npy2tsv")
    a.add_argument("source", help="Source .npy file")
    a.add_argument("destination", help="Location to save matrix as a .tsv file")
    args = a.parse_args()
    main(args.source, args.destination)
