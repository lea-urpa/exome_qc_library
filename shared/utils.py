"""
utility functions
"""
import os
from functools import partial
from tempfile import NamedTemporaryFile
import gzip
import csv
import subprocess
import shlex
import time
import logging
import hail as hl


def check_exists(filename):
    """
    Check if named file exists, adding needed slash if .ht or .mt object
    :param filename:
    :return:
    """
    if filename.startswith("gs://"):
        if filename.endswith(".mt") or filename.endswith(".ht"):
            filename = filename + "/"
        if filename.endswith(".mt/") or filename.endswith(".ht/"):
            filename = filename + "metadata.json.gz" # needs file, not mt/ht dir
        qstat_cmd = ['gsutil', '-q', 'stat', filename]
        exists = subprocess.call(qstat_cmd)

        if exists == 0:
            return True
        else:
            return False

    else:
        return os.path.exists(filename)


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


def copy_logs_output(log_dir, log_file, plot_dir):
    if not log_dir.endswith("/"):
        log_dir = log_dir + "/"

    datestr = time.strftime("%Y.%m.%d")
    hail_log_name = os.path.join(log_dir, datestr + "_hail_log.txt")
    hl.copy_log(hail_log_name)

    cmd = ['gsutil', 'cp', log_file, log_dir]
    subprocess.call(cmd)
    cmd = ['gsutil', 'cp', '*.html', plot_dir]
    subprocess.call(cmd)
    cmd = ['gsutil', 'cp', '*.pdf', plot_dir]
    subprocess.call(cmd)


def add_secondary(cluster_name, num_preemptibles, region='europe-west1'):
    """
    Add preemptible nodes on the cluster for applicable (non-shuffle) steps.
    :param cluster_name: Name of the cluster
    :param num_preemptibles: Number of preemptibles
    :param region: region (default europe-west1)
    :return:
    """
    if cluster_name is not None:
        logging.info(f"Adding {str(num_preemptibles)} preemptiple nodes to cluster {cluster_name}")

        cmd = shlex.split(f"gcloud dataproc clusters update {cluster_name} --region {region} "
                          f"--num-secondary-workers {str(num_preemptibles)}")

        subprocess.call(cmd)


def remove_secondary(cluster_name, region='europe-west1'):
    """
    Removes all preemptible nodes from a cluster.
    :param cluster_name: cluster name
    :param region: region, default 'europe-west1'
    :return:
    """
    if cluster_name is not None:
        logging.info('Removing preemptible nodes.')

        cmd = shlex.split(f"gcloud dataproc clusters update {cluster_name} --region {region} "
                          f"--num-secondary-workers 0")

        subprocess.call(cmd)


def check_counts(mt, args):
    counts = mt.count()
    defined_rows = mt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt.locus)))
    defined_cols = mt.aggregate_cols(hl.agg.count_where(hl.is_defined(mt.s)))
    if (counts[0] != args.start_rows) or (counts[1] != args.start_cols):
        print(f"Error! counts mismatch! Start rows: {args.start_rows}, start cols: {args.start_cols}, count: {counts}")
    else:
        print("counts OK")

    if (defined_rows != args.start_rows) or (defined_cols != args.start_cols):
        print(f"Error! defined mismatch! Start rows: {args.start_rows}, start cols: {args.start_cols}, "
              f"defined rows: {defined_rows}, defined cols: {defined_cols}")
    else:
        print("defined counts OK ")


def tmp_bash(cmd, check=False):
    """
    Creates a temporary bash file to run subprocess command files, useful when you want to have the output piped to
    another file
    :param cmd: command to run
    :param check:
    :return:
    """
    script_file = NamedTemporaryFile(delete=True)

    with open(script_file.name, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(cmd + "\n")
    os.chmod(script_file.name, 0o777)
    script_file.file.close()

    if check:
        subprocess.check_call(script_file.name)
    else:
        subprocess.call(script_file.name, stderr=subprocess.DEVNULL)
