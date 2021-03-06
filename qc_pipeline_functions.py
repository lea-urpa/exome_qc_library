"""
Wrapper functions for running exome QC pipeline with Hail.

Author: Lea Urpa, August 2020
"""
import time
import sys
import os
import logging
import hail as hl
from bokeh.io import output_file, save

# Import scripts
scripts_dir = sys.argv[sys.argv.index('--scripts_dir') + 1]  # Finds the scripts_dir arg and grabs the next string
scripts = ["variant_qc.py", "samples_qc.py", "samples_annotation.py", "variant_annotation.py", "helper_scripts.py"]
for script in scripts:
    hl.spark_context().addPyFile(os.path.join(scripts_dir, script))

import samples_annotation as sa
import helper_scripts as h
import variant_annotation as va
import variant_qc as vq
import samples_qc as sq


def load_data(args):
    '''
    Loads in dataset from given directory, either test or full datset, or creates new dataset.

    :param args: command line arguments
    :return: Returns matrix table if args.checkpoint == cpcounter (1), or else None
    '''
    if args.checkpoint != args.cpcounter:
        mt = None
        return mt

    h.add_preemptibles(args.cluster_name, args.num_preemptible_workers)

    datestr = time.strftime("%Y.%m.%d")
    mt = hl.read_matrix_table(args.mt)
    mt.annotate_globals(original_mt_input={'file': args.mt, 'date': datestr})

    try:
        test = hl.is_defined(mt.vep)
    except Exception as e:
        logging.error("Error! Input matrix table has not been VEP annotated!")
        logging.error(e)
        exit()

    if args.test:
        logging.info('Test flag given, filtering to on chrom 22.')
        if args.reference_genome == "GRCh38":
            chrom_codes = hl.array(["chr22", "chrX"])
        else:
            chrom_codes = hl.array(["22", "X"])

        mt = mt.filter_rows(chrom_codes.contains(mt.locus.contig))

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
        checkpoint_name = f"{step_info['prefix']}_{args.out_name}_{step_info['suffix']}_test.mt/"
    else:
        checkpoint_name = f"{step_info['prefix']}_{args.out_name}_{step_info['suffix']}.mt/"

    logging.info(f"Writing checkpoint {args.cpcounter}: {step}")

    step_name = f"{args.cpcounter}_{step}"
    mt = mt.annotate_globals(**{step_name: datestr})
    mt = mt.checkpoint(os.path.join(args.out_dir, checkpoint_name), overwrite=True)

    logging.info(f"Copying logs.")
    h.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

    return mt


def load_checkpoint(checkpoint, step, args):
    """
    Loads in data file from a given checkpoint
    :param checkpoint: checkpoint to load data from
    :param step: name of the step
    :param args: arguments namespace object
    :return: returns loaded matrix table
    """
    step_info = getattr(args, step)
    if args.test:
        checkpoint_name = f"{step_info['prefix']}_{args.out_name}_{step_info['suffix']}_test.mt/"
    else:
        checkpoint_name = f"{step_info['prefix']}_{args.out_name}_{step_info['suffix']}.mt/"

    h.add_preemptibles(args.cluster_name, args.num_preemptible_workers)

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
    if (args.checkpoint > args.cpcounter) | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        return mt

    # Data has been loaded with load_data above, this is the first step in the QC pipeline.
    logging.info('Annotating samples.')
    step = 'samples_annotation'

    ###################################################
    # Annotate with optional samples annotation files #
    ###################################################
    if args.samples_annotation_files is not None:
        annotation_files = args.samples_annotation_files.strip().split(",")

        for file in annotation_files:
            mt = sa.annotate_cols_from_file(mt, file, args.samples_delim, args.samples_col, args.samples_miss)

    ##############################
    # Annotate with bam metadata #
    ##############################
    if args.bam_metadata is not None:
        mt = sa.annotate_cols_from_file(mt, args.bam_metadata, args.bam_delim, args.bam_sample_col, args.bam_miss)

    ##############################################
    # Check that particular needed columns exist #
    ##############################################
    columns = ['chimeras_col', 'contamination_col']

    if args.fam_id is not None:
        columns.append('fam_id')
    if args.pat_id is not None:
        columns.append('pat_id')
    if args.mat_id is not None:
        columns.append('mat_id')
    if args.batch_col_name is not None:
        columns.append('batch_col_name')
    if args.pheno_col is not None:
        columns.append('pheno_col')

    for colname in columns:
        col = getattr(args, colname)
        try:
            test = hl.is_defined(mt[col])
        except Exception as e:
            logging.error(f"Error! Given column annotation {col} does not actually exist after inputting sample "
                          f"annotations.")
            logging.error(e)
            exit(1)

    if args.overwrite_checkpoints:
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
    if (args.checkpoint > args.cpcounter) | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        return mt

    step = "sample_removal"
    # Load data from after sample annotation, if we are starting at this checkpoint, else pass from prev step
    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'samples_annotation', args)

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

    if (args.checkpoint > args.cpcounter) | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        return mt

    step = "low_pass_variant_qc"

    # Load data from after sample removal, if we are starting at this checkpoint, else pass from prev step
    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'sample_removal', args)

    ########################################################
    # Annotate variants and genotypes for those failing QC #
    ########################################################
    logging.info("Running low-pass variant QC and genotype QC before samples QC.")
    mt = vq.find_failing_variants(mt, args, mode='low_pass')

    if args.overwrite_checkpoints:
        mt = save_checkpoint(mt, step, args)

    args.cpcounter += 1
    return mt


def maf_prune_relatedness(mt, args):
    """
    MAF prunes dataset for relatedness calculations in King.

    :param mt: matrix table to prune
    :param args:
    :return:
    """
    if (args.checkpoint > args.cpcounter) | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        mt_ldpruned = None
        return mt, mt_ldpruned

    step = "maf_ld_prune"

    # Load data from after phenotype samples QC, if we are starting at this checkpoint, else pass from prev step
    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'low_pass_variant_qc', args)

    #######################################################
    # Filter out failing samples, variants, and genotypes #
    #######################################################
    mt_filtered = sq.filter_failing(mt, args, mode='low_pass', unfilter_entries=True, checkpoint=False)

    ###############################
    # Filter out low MAF variants #
    ###############################
    mt_maffilt = vq.maf_filter(mt_filtered, args.ind_maf, filter_ac0_after_pruning=True)

    ##########################################
    # Downsample variants if row count > 80k #
    ##########################################
    mt_ldpruned = vq.downsample_variants(mt_maffilt, 80000)

    if args.overwrite_checkpoints:
        mt_ldpruned = save_checkpoint(mt_ldpruned, step, args)

    args.cpcounter += 1
    return mt, mt_ldpruned


def find_related_individuals(mt, mt_ldpruned, args):
    """
    Either exports genotype data as Plink files for King relatedness calculation, or annotates the mt based on
    previously calculated relatedness exclusions, depending on param.run_king.

    :param mt: matrix table to annotate with relatedness exclusion info
    :param mt_ldpruned: matrix table to export to Plink to run King
    :return: returns mt, mt_ldpruned
    """

    if (args.checkpoint > args.cpcounter) | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        return mt, mt_ldpruned

    step = "find_related_individuals"

    # Load data from after maf pruning, if we are starting at this checkpoint, else pass from prev step
    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'low_pass_variant_qc', args)
        mt_ldpruned = load_checkpoint(args.checkpoint, 'maf_ld_prune', args)

    ####################################
    # Export data as Plink to run King #
    ####################################
    if args.run_king is True:
        logging.info("Exporting MAF filtered and LD pruned dataset to Plink, for running King.")
        king_dir = os.path.join(args.out_dir, "king")

        outputs = {'ind_id': mt_ldpruned.s}
        if args.pheno_col is not None:
            outputs['pheno'] = mt_ldpruned[args.pheno_col]
        if args.fam_id is not None:
            outputs['fam_id'] = mt_ldpruned[args.fam_id]
        if args.pat_id is not None:
            outputs['pat_id'] = mt_ldpruned[args.pat_id]
        if args.mat_id is not None:
            outputs['mat_id'] = mt_ldpruned[args.mat_id]

        hl.export_plink(mt_ldpruned,  os.path.join(king_dir, args.out_name), **outputs)

        logging.info(f"Plink dataset exported, time to run King now and come back. Start agian at checkpoint:"
                     f"{args.checkpoint}")

        args.cpcounter += 1
        return mt, mt_ldpruned

    ######################################
    # Or load in relatedness information #
    ######################################
    else:
        try:  # Try annotating relateds, if it fails exit
            logging.info('Uploading list of relatives to remove + annotating matrix table.')
            logging.info('Note: not removing these individuals from the dataset, but marking them for removal in '
                         'steps such as PCA calculation. They are retained in the dataset.')
            mt = sa.annotate_relateds(mt, args.relatives_removal_file)
            mt_ldpruned = sa.annotate_relateds(mt_ldpruned, args.relatives_removal_file)

            if args.overwrite_checkpoints:
                mt = save_checkpoint(mt, step, args)
                mt_ldpruned = save_checkpoint(mt_ldpruned, step + "_ld_pruned", args)

            args.cpcounter += 1
            return mt, mt_ldpruned

        except Exception as e:
            logging.error('Annotating related individuals failed. Have you run King? Exiting now.')
            logging.info(e)

            args.run_king = True

            args.cpcounter += 1
            return mt, mt_ldpruned


def find_pop_outliers(mt, mt_ldpruned, args):
    """
    Find and annotate population outliers in the samples.

    :param mt: matrix table to annotate
    :param mt_ldpruned: LD pruned matrix table for calculating principal components
    :param args:
    :return: returns matrix table annotated with population outliers
    """

    if (args.checkpoint > args.cpcounter) | args.run_king | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        return mt, mt_ldpruned

    step = "find_pop_outliers"

    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'find_related_individuals', args)
        mt_ldpruned = load_checkpoint(args.checkpoint, 'find_related_individuals', args)

    ##########################################################################################
    # Find population outliers from filtered + LD pruned dataset, annotate unfiltered datset #
    ##########################################################################################
    h.remove_preemptibles(args.cluster_name)
    mt = sq.find_pop_outliers(mt_ldpruned, mt_to_annotate=mt, args=args)

    if args.overwrite_checkpoints:
        mt = save_checkpoint(mt, step, args)

    h.add_preemptibles(args.cluster_name, args.num_preemptible_workers)
    args.cpcounter += 1
    return mt


def impute_sex(mt, args):
    """
    Impute sex before variant QC filtering, run sex-specific annotations.

    :param mt:
    :param args:
    :return:
    """
    if (args.checkpoint > args.cpcounter) | args.run_king | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        return mt

    step = 'impute_sex'

    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'find_pop_outliers', args)

    ##########################################################
    # Run custom variant annotation, parsing VEP annotations #
    ##########################################################
    mt = va.annotate_variants(mt)

    #################################################################################
    # Filter out failing samples, variants, and genotypes, filter out rare variants #
    #################################################################################
    mt_maf = sq.filter_failing(mt, args, mode='low_pass', unfilter_entries=False)
    mt_maf = vq.maf_filter(mt_maf, 0.05)

    ##############
    # Impute sex #
    ##############
    mt_maf, imputed_sex, mt = sq.impute_sex_plot(mt_maf, mt_to_annotate=mt, args=args)

    ########################################
    # Annotate  with sex-aware annotations #
    ########################################
    logging.info("Annotating sex-aware sample annotations, using dataset with failing samples, variants, "
                 "and genotypes filtered out.")
    mt_filtered = sq.filter_failing(mt, args, mode='low_pass', unfilter_entries=False)
    mt = sa.sex_aware_sample_annotations(mt_filtered, mt_to_annotate=mt, args=args)
    mt = va.sex_aware_variant_annotations(mt_filtered, mt_to_annotate=mt, args=args)

    if args.overwrite_checkpoints:
        mt = save_checkpoint(mt, step,  args)
    args.cpcounter += 1

    return mt


def samples_qc(mt, args):
    """
    Performs samples quality control on the dataset
    :param mt: matrix table to analyze
    :param args:
    :return: returns matrix table annotated with samples failing analytical samples QC
    """
    if (args.checkpoint > args.cpcounter) | args.run_king | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        return mt

    step = "samples_qc"

    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'impute_sex', args)

    ######################################################
    # Filter failing variants + genotypes for samples QC #
    ######################################################
    logging.info("Removing variants and genotypes failing low pass samples QC, and population outlier samples.")
    mt_filtered = sq.filter_failing(mt, args, mode='low_pass', unfilter_entries=False)

    ##################
    # Run samples QC #
    ##################
    mt = sq.samples_qc(mt_filtered, mt_to_annotate=mt, args=args)

    if args.overwrite_checkpoints:
        mt = save_checkpoint(mt, step, args)

    args.cpcounter += 1
    return mt


def final_variant_qc(mt, args):
    """
    Performs final variant QC on the data

    :param mt: matrix table to annotate failing variants
    :param args:
    :return: returns matrix table with variants failing final QC annotated
    """
    if (args.checkpoint > args.cpcounter) | args.run_king | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        return mt

    step = "final_variant_qc"

    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'impute_sex', args)

    logging.info('Running variant QC filtering.')
    mt = vq.find_failing_variants(mt, args, mode='final')

    if args.overwrite_checkpoints:
        mt = save_checkpoint(mt, step, args)

    args.cpcounter += 1
    return mt


def find_failing_variants_by_pheno(mt, args):
    """
    Finds variants failing on per-phenotype information- allelic balance and call rate, in cases and controls separately
    :param mt: matrix table to annotate
    :param args:
    :return: matrix table with new row annotation for variants failing per-phenotype QC
    """
    if (args.checkpoint > args.cpcounter) | args.run_king | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        return mt

    step = "filter_variants_by_phenotype"

    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'final_variant_qc', args)

    if args.pheno_col is not None:
        logging.info('Finding variants failing by phenotype.')
        mt = vq.find_variants_failing_by_pheno(mt, args)

    else:
        logging.info("Phenotype column not given, skipping filtering variants by phenotype.")

    if args.overwrite_checkpoints:
        mt = save_checkpoint(mt, step, args)
    args.cpcounter += 1

    return mt


def calculate_final_pcs(mt, args):
    """
    Calculates final principal components, after filtering out failing samples and variants, and annotates unfiltered
    matrix table. Removes relatives and projects them back with PC loadings.

    :param mt: matrix table on which sample and variant QC has been run
    :param args:
    :return: matrix table annotate with principal components.
    """

    if (args.checkpoint > args.cpcounter) | args.run_king | (args.checkpoint >= args.stop_checkpoint):
        args.cpcounter += 1
        return mt

    step = 'final_pc_calculation'

    if args.checkpoint == args.cpcounter:
        mt = load_checkpoint(args.checkpoint, 'filter_variants_by_phenotype', args)
    #######################################################################
    # Filter out failing samples, variants, genotypes for PC calculations #
    #######################################################################
    mt_filtered = sq.filter_failing(mt, args, mode="final", unfilter_entries=True, pheno_qc=True)

    ########################
    # Filter out relatives #
    ########################
    logging.info("Filtering to unrelated individuals for PC calculations.")
    mt_norelateds = mt_filtered.filter_cols(mt_filtered.related_to_remove == False, keep=True)

    ##############################################
    # MAF prune and LD prune for calculating PCS #
    ##############################################
    logging.info("Filtering to common variants and LD pruning dataset.")
    mt_mafpruned = vq.maf_filter(mt_norelateds, 0.05)
    mt_mafpruned = save_checkpoint(mt_mafpruned, 'maf_filter_pcs', args)

    logging.info('LD pruning final dataset for PC calculation')
    mt_ldpruned = vq.ld_prune(mt_mafpruned, args)
    mt_ldpruned = save_checkpoint(mt_ldpruned, 'ld_prune_pcs', args)

    ################################################
    # Calculate PCs and project to relatives, plot #
    ################################################
    h.remove_preemptibles(args.cluster_name)
    mt = sq.project_pcs_relateds(mt_ldpruned, mt, args.pc_num)
    datestr = time.strftime("%Y.%m.%d")

    if args.pca_plot_annotations is not None:
        pca_annotations = args.pca_plot_annotations.strip().split(",")
        for annotation in pca_annotations:
            output_file(f"{datestr}_final_pcs_plot_{annotation}.html")
            p = hl.plot.scatter(mt.pc1, mt.pc2, label=mt[annotation])
            save(p)
    else:
        output_file(f'{datestr}_final_pcs_plot.html')
        pcplot = hl.plot.scatter(mt.pc1, mt.pc2)
        save(pcplot)

    if args.overwrite_checkpoints:
        mt = save_checkpoint(mt, step, args)

    args.cpcounter += 1
    return mt
