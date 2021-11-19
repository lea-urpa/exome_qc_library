"""
Short script for pulling Ensembl annotations with additional info in one column to separate columns-
parsing this is pretty slow in Hail, let's see if we can't speed it up by preprocessing the data here.

Author: Lea Urpa
Date: November 2021
"""
from ordered_set import OrderedSet
import os
import argparse

def expand_ensembl_gtf(filen):
    # Do first pass through to get all the headings from the annotations column
    colnames = OrderedSet(['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame'])

    print(f"Scanning file {filen} for headers.")
    with open(filen) as in_f:
        for line in in_f:
            if not line.startswith("#"):
                annotations = line.strip().split("\t")[8].split(";")
                annotations = [x.strip().split(" ") for x in annotations]
                headers = [x[0] for x in annotations if x[0] is not '']

                [colnames.add(x) for x in headers]
    in_f.close()

    # Pass through again and expand annotations to separate columns
    new_fn = filen.replace(".gtf", "") + "_expanded.tsv"
    print(f"Expanding annotations column for file {filen}. Writing to new file{new_fn}")
    new_gtf = open(new_fn, "w")
    new_gtf.write("\t".join(list(colnames)) + "\n")

    with open(filen) as in_f:
        linecount = 0
        for line in in_f:
            linecount += 1
            if not line.startswith("#"):
                # Split by word, pull annotations, split again
                words = line.strip().split("\t")
                annotations = words[8].strip().split(";")
                annotations = [x.strip().split(" ") for x in annotations if x is not ""]
                del words[-1]

                # Extend columns with fixed size empty NAs
                words.extend(["NA"] * (len(colnames) - len(words)))

                for annot in annotations:
                    annot_name = annot[0]
                    annot_value = annot[1].replace('"', '')
                    words[colnames.index(annot_name)] = annot_value

                new_gtf.write("\t".join(words) + "\n")

    in_f.close()
    new_gtf.close()

if __name__ == "__main__":
    ###################
    # Parse arguments #
    ###################
    parser = argparse.ArgumentParser(description="Expand ensembl GTF files to sparse TSVs for reading into R.")
    parser.add_argument("--gtf_filename", type=str, help="Name and file path of GTF file to expand.", required=True)

    args = parser.parse_args()

    # Expand GTF
    if os.path.isfile(args.gtf_filename):
        expand_ensembl_gtf(args.gtf_filename)
    else:
        print(f"Error! File {args.gtf_filename} not found. Exiting.")
