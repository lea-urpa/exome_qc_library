"""
Functions for use in samples QC for exome sequencing data with Hail.

Author: Lea Urpa, August 2020
"""
import hail as hl
import logging


def maf_prune(mt, args, varqc_name, filter_after_pruning=False, unfilter_entries=False, ):
    """
    Takes matrix table, filters out failing genotypes, variants, and samples, and MAF prunes the
    table, and returns the matrix table

    :param mt: matrix table to prune (should be LD pruned and have x chrom removed).
    :param filter_after_pruning: filter variants no longer in the data, e.g. sum(AC) = 0? re-runs variant qc.
    :param unfilter_entries: Unfilter entries following maf pruning? Necessary to run PCA later.
    :return: returns maf-pruned matrix table.
    """
    logging.info("Removing failing + pop outlier samples, genotypes, and variants, "
                 "including genotypes failing on allelic balance.")
    start_count = mt.count()

    mt = mt.filter_entries((hl.len(mt.failing_depth_quality) == 0) & hl.len(mt.failing_ab) == 0, keep=True)
    mt = mt.filter_samples((hl.len(mt.failing_samples_qc) == 0) & mt.population_outlier == False, keep=True)
    mt = mt.filter_variants(hl.len(mt[varqc_name]) == 0, keep=True)

    final_count = mt.count()

    logging.info(f"Matrix table count before filtering: {start_count}. After filtering: {final_count}")

    # Run hl.variant_qc() to get AFs
    mt = hl.variant_qc(mt)

    # Filter MAF
    logging.info(f'Filtering out variants with minor allele frequency < {args.ind_maf}')
    mt = mt.filter_rows(mt.row.variant_qc.AF[1] > args.ind_maf, keep=True)
    mt = mt.annotate_globals(maf_threshold_LDpruning=args.ind_maf)

    if filter_after_pruning:
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(hl.sum(mt.row.variant_qc.AC) == hl.int(0), keep=False)

    if unfilter_entries is True:  # Unfilter entries if needed for pc_relate
        mt = mt.unfilter_entries()
        logging.info('Unfiltered entries (set to missing) successfully.')

    logging.info("MAF pruned mt count:" + str(mt.count()))
    return mt
