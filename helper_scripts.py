"""
Helper scripts for running Hail pipelines.

Author: Lea Urpa, August 2020
"""
import os
import time
import subprocess
import logging
import shlex
import hail as hl


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


def add_preemptibles(cluster_name, num_preemptibles, region='europe-west1'):
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


def remove_preemptibles(cluster_name, region='europe-west1'):
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
