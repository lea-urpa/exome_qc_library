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
import samples_qc as sq


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

        mt = mt.filter_rows(mt.locus.contig.contains(chrom_codes))

        checkpoint_name = args.mt.replace(".mt", "") + "_test.mt"
        mt = mt.checkpoint(checkpoint_name, overwrite=True)

    return mt


def save_checkpoint(mt, step, args, pruned=False):
    '''
    Saves a matrix table at a checkpoint with custom global annotation and file ending.
    :param mt: matrix table to checkpoint
    :param step: step the pipeline is on, e.g. low pass variant QC
    :param args: arguments namespace object
    :return: returns checkpointed matrix table
    '''
    datestr = time.strftime("%Y.%m.%d")
    step_info = getattr(args, step)

    if pruned:
        prune_str = "_filtered_pruned"
    else:
        prune_str = ""

    if args.test:
        checkpoint_name = f"{step_info['prefix']}_{args.out_name}_{step_info['suffix']}{prune_str}_test.mt"
    else:
        checkpoint_name = f"{step_info['prefix']}_{args.out_name}_{step_info['suffix']}{prune_str}_.mt"

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


def remove_samples(mt, args):
    """
    Removes samples based on arbitrary lists and key words.

    :param mt:
    :param args:
    :return:
    """
    if args.checkpoint > args.cpcounter:
        args.cpcounter += 1
        return mt

    step = "sample_removal"
    # Load data from after sample annotation, if we are starting at this checkpoint, else pass from prev step
    if args.checkpoint == args.cpcounter:
        mt = hl.load_checkpoint(args.checkpoint, 'samples_annotation', args)

    if (args.sample_removal_strings is not None) or (args.sample_removal_list is not None):
        samples_start = mt.count_cols()
        logging.info(f"Initial sample count: {samples_start}")

    # Filter out samples from arbitrary list uploaded with args
    if args.sample_removal_list is not None:
        rm_list = hl.import_table(args.sample_removal_list, no_header=True)
        rm_list = rm_list.annotate(s=rm_list.f0)
        rm_list = rm_list.key_by("s")
        mt = mt.anti_join_cols(rm_list)
        list_filtered = mt.count_cols()
        logging.info(f"Sample count after filtering by input list: {list_filtered}")

    # Filter out samples that have sample names containing a particular string
    if args.sample_removal_strings is not None:
        removal_strings = args.sample_removal_strings.strip().split(",")
        for key_string in removal_strings:
            mt = mt.filter_cols(mt.s.contains(key_string), keep=False)
        string_filtered = mt.count_cols()
        logging.info(f"Sample count after filtering by strings: {string_filtered}")

    if args.overwrite_checkpoints:
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

    # Load data from after sample removal, if we are starting at this checkpoint, else pass from prev step
    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'sample_removal', args)

    # Add preemtible nodes
    h.add_preemptibles(args.cluster_name, args.num_preemptible_workers)

    ########################################################
    # Annotate variants and genotypes for those failing QC #
    ########################################################
    logging.info("Running low-pass variant QC and genotype QC before samples QC.")
    mt = vq.find_failing_variants(mt, args, mode='low_pass')

    # Remove the preemptible nodes
    h.remove_preemptibles(args.cluster_name)

    if args.overwrite_checkpoints:
        mt = save_checkpoint(mt, step, args)

    args.cpcounter += 1
    return mt


def maf_LDprune_relatedness(mt, args):
    """
    MAF prunes dataset for relatedness calculations in King.

    :param mt: matrix table to prune
    :param args:
    :return:
    """
    if args.checkpoint > args.cpcounter:
        args.cpcounter += 1
        mt_mafpruned = None
        return mt, mt_mafpruned

    step = "maf_prune_relatedness"

    # Load data from after phenotype samples QC, if we are starting at this checkpoint, else pass from prev step
    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'low_pass_variant_qc', args)

    #######################################################
    # Filter out failing samples, variants, and genotypes #
    #######################################################
    mt_filtered = sq.filter_failing(mt, args, varqc_name=args.lowpass_fail_name, unfilter_entries=False)

    h.add_preemptibles(args.cluster_name, args.num_preemptible_workers)

    ###############################
    # Filter out low MAF variants #
    ###############################
    mt_maffilt = sq.maf_filter(mt_filtered, args,  filter_ac0_after_pruning=True)

    ################################################
    # LD prune if there are more than 80k variants #
    ################################################
    if not mt_maffilt.row_count() < 80000:
        mt_ldpruned = sq.ld_prune(mt_maffilt, args)
    else:
        logging.info("Detected MAF filtered dataset already has less than 80 000 variants, skipping LD pruning for King"
                     "relatedness calculations.")
        mt_ldpruned = mt_maffilt

    h.remove_preemptibles(args.cluster_name)

    if args.overwrite_checkpoints:
        mt_ldpruned = save_checkpoint(mt_ldpruned, step, args, pruned=True)
    args.cpcounter += 1
    return mt, mt_ldpruned


def find_related_individuals(mt, mt_mafpruned, args):
    """
    Either exports genotype data as Plink files for King relatedness calculation, or annotates the mt based on
    previously calculated relatedness exclusions, depending on param.run_king.

    :param mt: matrix table to annotate with relatedness exclusion info
    :param mt_mafpruned: matrix table to export to Plink to run King
    :return: returns mt, mt_mafpruned
    """

    if args.checkpoint > args.cpcounter:
        args.cpcounter += 1
        return mt, mt_mafpruned

    step = "find_related_individuals"

    # Load data from after maf pruning, if we are starting at this checkpoint, else pass from prev step
    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'low_pass_variant_qc', args)
        mt_mafpruned = load_checkpoint(args.checkpoint, 'maf_prune_relatedness', args)

    ####################################
    # Export data as Plink to run King #
    ####################################
    if args.run_king is True:
        logging.info("Exporting MAF filtered and LD pruned dataset to Plink, for running King.")
        king_dir = os.path.join(args.out_dir, "king")

        outputs = {'ind_id': mt_mafpruned.s}
        if args.pheno_col is not None:
            outputs['pheno'] = mt_mafpruned[args.pheno_col]
        if args.fam_id is not None:
            outputs['fam_id'] = mt_mafpruned[args.fam_id]
        if args.pat_id is not None:
            outputs['pat_id'] = mt_mafpruned[args.pat_id]
        if args.mat_id is not None:
            outputs['mat_id'] = mt_mafpruned[args.mat_id]

        hl.export_plink(mt_mafpruned,  os.path.join(king_dir, args.out_name), **outputs)

        logging.info(f"Plink dataset exported, time to run King now and come back. Start agian at checkpoint:"
                     f"{args.checkpoint}")

        args.cpcounter += 1
        return mt, mt_mafpruned

    ######################################
    # Or load in relatedness information #
    ######################################
    else:
        try:  # Try annotating relateds, if it fails exit
            logging.info('Uploading King relatedness information + annotating matrix table.')
            mt = sa.annotate_relateds(mt, args.relatives_file)
            mt_mafpruned = sa.annotate_relateds(mt_mafpruned, args.relatives_file)
            if args.overwrite_checkpoints:
                mt = save_checkpoint(mt, step, args)
                mt_mafpruned = save_checkpoint(mt_mafpruned, step, args, pruned=True)

            args.cpcounter += 1
            return mt, mt_mafpruned

        except Exception as e:
            logging.error('Annotating related individuals failed. Have you run King? Exiting now.')
            logging.info(e)

            args.run_king = True

            args.cpcounter += 1
            return mt, mt_mafpruned