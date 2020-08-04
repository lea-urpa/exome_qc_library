"""
Helper scripts for running Hail pipelines.

Author: Lea Urpa, August 2020
"""
import hail as hl
import logging
import time
import sys
import subprocess


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
