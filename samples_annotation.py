"""
Functions for annotating samples in a matrix table with Hail.

Author: Lea Urpa, August 2020
"""
import hail as hl


def annotate_cols_from_file(mt, file, args):

    if not file.endswith(".ht"):
        ht = hl.import_table(file, delimiter=args.samples_delim, impute=True, key=args.samples_col,
                             missing=args.samples_miss)
    else:
        ht = hl.read_table(file)


    # TODO next: figure out how to take all the columns of a hail table, minus the key, and use that
    # to annotate the columns of the matrix table