"""
Helper scripts for running Hail pipelines.

Author: Lea Urpa, August 2020
"""
import time
import subprocess
import logging
import shlex
import hail as hl


def configure_logging(logstem):
    '''
    Configures logging for pipelines submitted via `cluster submit mycluster myscript.py`

    :param logstem: string giving some reasonable description of the pipeline's purpose, for use in log file names
    :return: Returns datestring, timestring, and string with log file name.
    '''
    # Define time and date strings
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day

    log_file = logstem + '_' + timestr + '.txt'

    return datestr, timestr, log_file


def copy_logs_output(log_dir, timestr, log_file, plot_dir):
    hl.copy_log(log_dir + "hail_log" + timestr + ".txt")
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
