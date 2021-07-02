"""
Wrapper functions for running exome QC pipeline with Hail.

Author: Lea Urpa, August 2020
"""
import time
import sys
import os
import logging
import hail as hl
import samples_annotation as sa




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



