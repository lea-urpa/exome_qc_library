"""
utility functions
"""
import os
from functools import partial
import gzip
import csv

def get_path_info(path):  # shamelessly stolen from Pietro
    file_path = os.path.dirname(path)
    basename = os.path.basename(path)
    file_root, file_extension = os.path.splitext(basename)
    return file_path, file_root, file_extension


def return_open_func(f):  # shamelessly stolen from Pietro
    """
    Detects file extension and return proper open_func
    """

    file_path, file_root, file_extension = get_path_info(f)

    if 'bgz' in file_extension:
        open_func = partial(gzip.open, mode='rb')
    elif 'gz' in file_extension:
        open_func = partial(gzip.open, mode='rt')
    else:
        open_func = open
    return open_func


def check_delimiter(f):  # shamelessly stolen from Pietro
    open_func = return_open_func(f)
    with open_func(f) as i:header = i.readline().strip()
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(header)
    if dialect.delimiter != "0":
        return dialect.delimiter
    else:
        return None
