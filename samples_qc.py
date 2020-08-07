"""
Functions for use in samples QC for exome sequencing data with Hail.

Author: Lea Urpa, August 2020
"""
import hail as hl
import logging


def filter_failing(mt, args, varqc_name, unfilter_entries=False):
    """

    :param mt:
    :param args:
    :param varqc_name:
    :param unfilter_entries: Unfilter entries following maf pruning? Necessary to run PCA later.
    :return:
    """
    logging.info("Removing failing + pop outlier samples, genotypes, and variants, "
                 "including genotypes failing on allelic balance.")
    start_count = mt.count()

    mt = mt.filter_entries((hl.len(mt.failing_depth_quality) == 0) & hl.len(mt.failing_ab) == 0, keep=True)
    mt = mt.filter_samples((hl.len(mt.failing_samples_qc) == 0) & mt.population_outlier == False, keep=True)
    mt = mt.filter_variants(hl.len(mt[varqc_name]) == 0, keep=True)

    final_count = mt.count()

    logging.info(f"Matrix table count before filtering: {start_count}. After filtering: {final_count}")

    if unfilter_entries is True:  # Unfilter entries if needed for pc_relate
        logging.info('Unfiltering entries (setting them to missing)')
        mt = mt.unfilter_entries()

    logging.info(f"Writing temporary checkpoint for filtered mt.")
    mt = mt.checkpoint(f"{args.outputstem}_{varqc_name}_removed_tmp.mt")

    return mt


def maf_filter(mt, args, filter_ac0_after_pruning=False):
    """
    Takes matrix table, filters out failing genotypes, variants, and samples, and MAF prunes the
    table, and returns the matrix table

    :param mt: matrix table to prune (should be LD pruned and have x chrom removed).
    :param filter_ac0after_pruning: filter variants no longer in the data, e.g. sum(AC) = 0?
    :return: returns maf filtered matrix table.
    """
    # Run hl.variant_qc() to get AFs
    mt = hl.variant_qc(mt)

    # Filter MAF
    logging.info(f'Filtering out variants with minor allele frequency < {args.ind_maf}')
    mt = mt.filter_rows(mt.row.variant_qc.AF[1] > args.ind_maf, keep=True)
    mt = mt.annotate_globals(maf_threshold_LDpruning=args.ind_maf)

    if filter_ac0_after_pruning:
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(hl.sum(mt.row.variant_qc.AC) == hl.int(0), keep=False)

    logging.info("MAF pruned mt count:" + str(mt.count()))
    return mt


def ld_prune(mt, args, rm_chr_x=False):
    '''
    LD prune and remove chromosome X from a matrix table, for calculating kinship and principal components

    :param mt: matrix table to annotate, should already have related individuals removed.
    :return: returns the ld pruned matrix table
    '''

    # LD prune
    pruned_variant_table = hl.ld_prune(mt.GT, r2=args.r2, bp_window_size=args.bp_window_size)
    mt_ldpruned = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))

    # Remove chromosome X
    if rm_chr_x:
        if args.reference_genome == "GRCh38":
            chrom = "chrX"
        elif args.reference_genome == "GRCh37":
            chrom = "X"
        mt_ldpruned = hl.filter_rows(mt.locus.contig == chrom, keep=False)

    logging.info(f"Variant and sample count after LD pruning: {mt_ldpruned.count()}")

    mt_ldpruned = mt_ldpruned.annotate_globals(ld_pruning_parameters={'r2': args.r2,
                                                                      'bp_window_size': args.bp_window_size})

    return mt_ldpruned
