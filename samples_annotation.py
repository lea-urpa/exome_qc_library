"""
Functions for annotating samples in a matrix table with Hail.

Author: Lea Urpa, August 2020
"""
import hail as hl


def annotate_cols_from_file(mt, file, args):
    """
    Takes a matrix table, a file with sample information, and annotates the cols of the matrix table with each column
    in the given file. File can be hail table, as long as the file ends with .ht it is detected and read.

    :param mt: matrix table with columns keyed with variable 's'
    :param file: text file to import with hl.import_table, or hail table to read with hl.read_table.
    :param args: arguments to get samples files delimiter, column with sample IDs, and missing data value.
    :return:
    """
    if not file.endswith(".ht"):
        ht = hl.import_table(file, delimiter=args.samples_delim, impute=True, key=args.samples_col,
                             missing=args.samples_miss)
    else:
        ht = hl.read_table(file)

    cols = list(dict(ht.row).keys())
    cols.remove(args.samples_col)

    for colname in cols:
        mt = mt.annotate_cols(**{colname: ht[mt.s][colname]})

    return mt
