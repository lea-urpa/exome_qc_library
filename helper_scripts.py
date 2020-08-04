"""
Helper scripts for running Hail pipelines.

Author: Lea Urpa, August 2020
"""
import hail as hl
import logging
import time
import sys
import subprocess


def configure_logging(logstem, logging_level=logging.INFO):
    '''
    Configures logging for pipelines submitted via `cluster submit mycluster myscript.py`

    :param logstem: string giving some reasonable description of the pipeline's purpose, for use in log file names
    :param logging_level: logging level- info by default
    :return: Returns datestring, timestring, and string with log file name.
    '''
    # Define time and date strings
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day

    log_file = logstem + '_' + timestr + '.txt'

    # Configure logger
    root = logging.getLogger()
    log_formatter = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(filename=log_file, format=log_formatter, level=logging_level)

    handler = logging.StreamHandler(sys.stdout)
    root.addHandler(handler)

    return datestr, timestr, log_file


def copy_logs_output(log_dir, timestr, log_file, out_dir):
    hl.copy_log(log_dir + "hail_log" + timestr + ".txt")
    cmd = ['gsutil', 'cp', log_file, log_dir]
    subprocess.call(cmd)
    cmd = ['gsutil', 'cp', '*.html', out_dir]
    subprocess.call(cmd)
    cmd = ['gsutil', 'cp', '*.pdf', out_dir]
    subprocess.call(cmd)
