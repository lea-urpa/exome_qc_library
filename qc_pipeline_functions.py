"""
Wrapper functions for running exome QC pipeline with Hail.

Author: Lea Urpa, August 2020
"""
import time
import sys
import os
import logging
import hail as hl

# Import scripts
scripts_dir = sys.argv[sys.argv.index('--scripts_dir') + 1]  # Finds the scripts_dir arg and grabs the next string
scripts = ["variant_qc.py", "samples_qc.py", "samples_annotation.py", "variant_annotation.py", "helper_scripts.py"]
for script in scripts:
    hl.spark_context().addPyFile(scripts_dir + script)

import samples_annotation as sa
import helper_scripts as h
import variant_qc as vq


def load_data(args):
    '''
    Loads in dataset from given directory, either test or full datset, or creates new dataset.

    :param args: command line arguments
    :return: Returns matrix table if args.checkpoint == cpcounter (0), or else None
    '''
    if args.checkpoint != args.cpcounter:
        mt = None
        return mt

    datestr = time.strftime("%Y.%m.%d")

    mt = hl.read_matrix_table(args.mt)
    mt.annotate_globals(original_mt_input={'file': args.mt, 'date': datestr})

    if args.test:
        logging.info('Test flag given, filtering to on chrom 22.')
        if args.reference_genome == "GRCh38":
            chrom_codes = ["chr22", "chrX"]
        else:
            chrom_codes = ["22", "X"]

        mt = mt.filter_rows(mt.locus.contig in chrom_codes)

        checkpoint_name = args.mt.replace(".mt", "") + "_test.mt"
        mt = mt.checkpoint(checkpoint_name, overwrite=True)

    return mt


def save_checkpoint(mt, step, args):
    '''
    Saves a matrix table at a checkpoint with custom global annotation and file ending.
    :param mt: matrix table to checkpoint
    :param step: step the pipeline is on, e.g. low pass variant QC
    :param args: arguments namespace object
    :return: returns checkpointed matrix table
    '''
    datestr = time.strftime("%Y.%m.%d")
    step_info = getattr(args, step)
    if args.test:
        checkpoint_name = f"{step_info['prefix']}_{args.out_name}_{step_info['suffix']}_test.mt"
    else:
        checkpoint_name = f"{step_info['prefix']}_{args.out_name}_{step_info['suffix']}_.mt"

    logging.info(f"Writing checkpoint after {args.cpcounter}: {step}")

    step_name = f"{args.cpcounter}_{step}"
    mt = mt.annotate_globals(**{step_name: datestr})
    mt = mt.checkpoint(os.path.join(args.out_dir, checkpoint_name), overwrite=True)

    return mt


def load_checkpoint(checkpoint, step, args):
    """
    Loads in data file from a given checkpoint
    :param checkpoint: checkpoint to load data from
    :param step: name of the step
    :param args: arguments namespace object
    :return: returns loaded matrix table
    """
    '''
    Loads a data file from a given checkpoint
    '''
    step_info = getattr(args, step)
    if args.test:
        checkpoint_name = f"{step_info['prefix']}_{args.out_name}_{step_info['suffix']}_test.mt"
    else:
        checkpoint_name = f"{step_info['prefix']}_{args.out_name}_{step_info['suffix']}.mt"

    mt = hl.read_matrix_table(os.path.join(args.out_dir, checkpoint_name))
    logging.info(f"Starting at checkpoint {str(checkpoint)}. Loading data from after {step}.")

    return mt


def annotate_samples(mt, args):
    """
    Annotate matrix table with various sample information. Data has been loaded with load_data function above, this is
    the first step in the QC pipeline.

    :param mt: matrix table to annotate
    :param args: arguments namespace object
    :return: returns matrix table + checkpoint counter
    """
    if args.checkpoint > args.cpcounter:
        args.cpcounter += 1
        return mt

    # Data has been loaded with load_data above, this is the first step in the QC pipeline.
    logging.info('Annotating samples.')
    step = 'samples_annotation'

    annotation_files = args.samples_annotation_files.strip().split(",")

    for file in annotation_files:
        mt = sa.annotate_cols_from_file(mt, file, args)

    if args.force:
        mt = save_checkpoint(mt, step, args)

    args.cpcounter += 1
    return mt


def low_pass_var_qc(mt, args):
    """
    Performs low pass variant QC on a matrix table (for before samples QC)
    :param mt: matrix table to filter
    :param args: arguments with threshold information
    :return: returns matrix table with bad variants filtered out.
    """

    if args.checkpoint > args.cpcounter:
        args.cpcounter += 1
        return mt

    step = "low_pass_variant_qc"

    # Load data from after sample annotation, if we are starting at this checkpoint, else pass from prev step
    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'samples_annotation', args)

    # Add preemtible nodes
    h.add_preemptibles(args.cluster_name, args.num_preemptible_workers)

    # Run Hail variant and samples QC to get metrics
    logging.info('Running variant and samples QC.')
    mt = hl.variant_qc(mt, name='prefilter_variant_qc')
    mt = hl.sample_qc(mt, name='prefilter_sample_qc')

    # Filter out bad variants and genotypes
    logging.info("Running low-pass variant QC and genotype QC before samples QC.")
    mt_filt = vq.lowpass_filter_variants(mt, args)

    # Remove the preemptible nodes
    h.remove_preemptibles(args.cluster_name)

    if args.overwrite_checkpoints:
        mt_filt = save_checkpoint(mt_filt, step, args)

    args.cpcounter += 1
    return mt_filt
