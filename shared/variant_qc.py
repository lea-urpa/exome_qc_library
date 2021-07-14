"""
Functions for variant QC for exome sequencing data with Hail.

Author: Lea Urpa, August 2020
"""
import logging
import hail as hl
import samples_qc as sq
import utils


def downsample_variants(mt, target_count):
    """
    Takes a matrix table, checks if the variant count is above a target count, and if so downsamples the matrix table.
    :param mt: matrix table to downsample
    :param target_count: target number of variants the matrix table should be
    :return: downsampled matrix table, if mt count > target count, or original matrix table if not.
    """
    var_count = mt.count_rows()
    if var_count > target_count:
        logging.info(f"Matrix table has more than {target_count} variants, randomly downsampling to {target_count} variants.")
        keep_fraction = target_count/var_count
        mt = mt.sample_rows(keep_fraction)

        count = mt.count()
        logging.info(f"Matrix table count after downsampling: {count}")

        return mt
    else:
        logging.info(f"Matrix table already has fewer than {target_count} variants, skipping downsampling.")

        return mt


def maf_filter(mt, maf, filter_ac0_after_pruning=False):
    """
    Takes matrix table, filters out failing genotypes, variants, and samples, and MAF prunes the
    table, and returns the matrix table

    :param mt: matrix table to prune (should be LD pruned and have x chrom removed).
    :param filter_ac0_after_pruning: filter variants no longer in the data, e.g. sum(AC) = 0?
    :return: returns maf filtered matrix table.
    """

    # Run hl.variant_qc() to get AFs
    mt = hl.variant_qc(mt)

    # Filter MAF
    logging.info(f'Filtering out variants with minor allele frequency < {maf}')
    mt = mt.filter_rows(mt.row.variant_qc.AF[1] > maf, keep=True)
    mt = mt.annotate_globals(maf_threshold_LDpruning=maf)

    if filter_ac0_after_pruning:
        logging.info('Removing variants with alt allele count = 0 (monomorphic variants).')
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(hl.sum(mt.row.variant_qc.AC) == hl.int(0), keep=False)
        count = mt.count()
        logging.info(f"MT count after removing monomorphic variants and MAF filtering: {count}")
    else:
        logging.info("MAF pruned mt count:" + str(mt.count()))

    return mt


def ld_prune(mt, args):
    """
     LD prune a matrix table, for calculating kinship and principal components

    :param mt: matrix table to annotate, should already have related individuals removed.
    :param args: namespace object with threshold arguments
    :return: returns the ld pruned matrix table
    """

    pruned_variant_table = hl.ld_prune(mt.GT, r2=args.r2, bp_window_size=args.bp_window_size)
    mt_ldpruned = mt.filter_rows(hl.is_defined(pruned_variant_table[mt.row_key]))

    logging.info(f"Variant and sample count after LD pruning: {mt_ldpruned.count()}")

    mt_ldpruned = mt_ldpruned.annotate_globals(
        ld_pruning_parameters={'r2': args.r2, 'bp_window_size': args.bp_window_size})

    return mt_ldpruned


def filter_failing_GTs_depth_quality(mt, checkpoint_name, prefix="", min_dp=10, min_gq=20, count_failing=True,
                                     filter_missing_measures=False):
    """
    Filter out genotypes with low depth and GQ values.

    :param mt: matrix table to filter
   :param checkpoint_name: Name of checkpoint to write to, including bucket name or file dir to write to
    :param prefix: prefix to add to variant annotations
    :param min_dp: minimum allowed sequencing depth
    :param min_gq: minimum allowed genotype quality
    :param count_failing: Count number of genotypes that are failing QC measures? Slow but useful for troubleshooting
    :param filter_missing_measures: Filter out genotypes that are missing QC measures, or keep them in?
    :return:
    """
    logging.info(f"Filtering genotypes with min_dp < {min_dp}, or min GQ < {min_gq}.")
    if (not prefix.endswith("_")) and (prefix != ""):
        prefix = prefix + "_"

    # Get starting genotypes count, instantiate genotype annotation
    gt_total = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    gt_het_homalt = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) &
                                                            (mt.GT.is_hom_var() | mt.GT.is_het())))
    gt_homref = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_hom_ref()))

    ############################
    # Define filter conditions #
    ############################
    failing_dp_cond = ((mt.DP < min_dp) & hl.is_defined(mt.DP))
    pl_cond = ((mt.GT.is_het() | mt.GT.is_hom_var()) & (mt.PL[0] < min_gq) &
                    hl.is_defined(mt.GT) & hl.is_defined(mt.PL))
    gq_cond = (mt.GT.is_hom_ref() & (mt.GQ < min_gq) & hl.is_defined(mt.GQ) & hl.is_defined(mt.GT))

    missing_dp_cond = (~(hl.is_defined(mt.DP)) & hl.is_defined(mt.GT))
    missing_gq_cond = (~(hl.is_defined(mt.GQ)) & hl.is_defined(mt.GT) & mt.GT.is_hom_ref())
    missing_pl_cond = (~(hl.is_defined(mt.PL)) & (mt.GT.is_het() | mt.GT.is_hom_var()) & hl.is_defined(mt.GT))

    ########################################
    # Count failing genotypes if indicated #
    ########################################
    if count_failing:
        # Count failing depth
        failing_depth = mt.aggregate_entries(hl.agg.count_where(failing_dp_cond))
        failing_depth_perc = round(failing_depth / gt_total*100, 2)
        missing_depth = mt.aggregate_entries(hl.agg.count_where(missing_dp_cond))
        missing_depth_perc = round(missing_depth / gt_total*100, 2)

        # Count failing GQ or PL measures
        failing_gq = mt.aggregate_entries(hl.agg.count_where(gq_cond))
        failing_pl = mt.aggregate_entries(hl.agg.count_where(pl_cond))
        missing_gq = mt.aggregate_entries(hl.agg.count_where(missing_gq_cond))
        missing_pl = mt.aggregate_entries(hl.agg.count_where(missing_pl_cond))

        logging.info(f"Number of GTs with depth < {min_dp}: {failing_depth} ({failing_depth_perc}%)")
        logging.info(f"Number of hom ref GTs with GQ < {min_gq}: {failing_gq} "
                     f"({round(failing_gq/gt_total*100, 2)}% of all GTs, "
                     f"{round(failing_gq/gt_homref*100, 2)}% of homref GTs)")
        logging.info(f"Number of het or hom alt GTs with PL[0] < {min_gq}: {failing_pl} "
                     f"({round(failing_pl/gt_total*100, 2)}% of all GTs, "
                     f"{round(failing_pl/gt_het_homalt*100, 2)} % of homalt + het GTs)")

        if missing_gq > 0:
            logging.info(f"Number of hom ref GTs that are defined but missing GQ : {missing_gq} "
                         f"({round(missing_gq/gt_total*100, 2)}% of all GTs, "
                         f"{round(missing_gq/gt_homref*100, 2)}% of homref GTs)")
            logging.info("(these genotypes are counted as failing)")
        if missing_pl > 0:
            logging.info(f"Number of het or hom alt GTs that are defined but missing PL: {missing_pl} "
                         f"({round(missing_pl/gt_total*100, 2)}% of all GTs, "
                         f"{round(missing_pl/gt_het_homalt*100, 2)}% of het and hom alt GTs)")
        if missing_depth > 0:
            logging.info(f"Number of GTs that are defined but missing depth measure: {missing_depth} ({missing_depth_perc}%)")
            logging.info("(these genotypes missing depth measure are counted as failing)")

        globals_fail_annot = prefix + "genotype_qc_failing_quality_depth"
        mt = mt.annotate_globals(
            **{globals_fail_annot: {'depth_excluded': failing_depth, 'quality_excluded': failing_gq}})

        #############################################################
        # Define filter conditions (again), or else source mismatch #
        #############################################################
        failing_dp_cond = ((mt.DP < min_dp) & hl.is_defined(mt.DP))
        pl_cond = ((mt.GT.is_het() | mt.GT.is_hom_var()) & (mt.PL[0] < min_gq) &
                   hl.is_defined(mt.GT) & hl.is_defined(mt.PL))
        gq_cond = (mt.GT.is_hom_ref() & (mt.GQ < min_gq) & hl.is_defined(mt.GQ) & hl.is_defined(mt.GT))

        missing_dp_cond = (~(hl.is_defined(mt.DP)) & hl.is_defined(mt.GT))
        missing_gq_cond = (~(hl.is_defined(mt.GQ)) & hl.is_defined(mt.GT) & mt.GT.is_hom_ref())
        missing_pl_cond = (~(hl.is_defined(mt.PL)) & (mt.GT.is_het() | mt.GT.is_hom_var()) & hl.is_defined(mt.GT))

    ############################
    # Filter failing genotypes #
    ############################
    if filter_missing_measures:
        logging.info("Filtering genotypes missing DP, GQ, or PL measures.")
        filter_condition = failing_dp_cond | pl_cond | gq_cond | missing_dp_cond |missing_gq_cond | missing_pl_cond
    else:
        logging.info("Filtering only genotypes failing QC measures, leaving genotypes missing those measures.")
        filter_condition = failing_dp_cond | pl_cond | gq_cond

    mt = mt.filter_entries(filter_condition)
    mt = mt.checkpoint(checkpoint_name + "_DP_GQ_filtered.mt/", overwrite=True)

    passing_gts = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    passing_gts_perc = round(passing_gts / gt_total * 100, 2)

    logging.info(f"Number of passing genotypes: {passing_gts} ({passing_gts_perc}%)")

    return mt


def count_variant_ab(mt, checkpoint_name, prefix="", samples_qc=False, pheno_col=None, min_het_ref_reads=0.2, max_het_ref_reads=0.8,
                     min_hom_ref_ref_reads=0.9, max_hom_alt_ref_reads=0.1):
    """
    :param mt: matrix table to annotate
    :param checkpoint_name: Name of checkpoint to write to, including bucket name or file dir to write to
    :param prefix: prefix to add to variant annotations
    :param samples_qc: has samples QC been run, failing samples and population outliers found?
    :param pheno_col: is there a column annotation giving True/False for case status?
    :param min_het_ref_reads: minimum percent reference reads for het genotype
    :param max_het_ref_reads: maximum percent reference reads for het genotype
    :param min_hom_ref_ref_reads: minimum percent reference reads for hom ref genotype
    :param max_hom_alt_ref_reads: maximum percent reference reads for hom alt genotype
    :return:
    """
    if (not prefix.endswith("_")) and (prefix != ""):
        prefix = prefix + "_"

    het_gt_cnt = prefix + 'n_het'
    het_ab = prefix + 'frac_het_gts_in_ab'
    case_het_gt_count = prefix + 'case_n_het'
    case_het_gt_ab = prefix + 'case_frac_het_gts_in_ab'
    cont_het_gt_count = prefix + 'control_n_het'
    cont_het_gt_ab = prefix + 'control_frac_het_gts_in_ab'

    ##########################################################################
    # Get total # of het GTs per variant in passing samples (or all samples) #
    ##########################################################################
    het_gt_filters = mt.GT.is_het() & hl.is_defined(mt.GT)

    if samples_qc:
        sample_filter = (
                (mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) == 0) &
                hl.is_defined(mt.pop_outlier_sample) & hl.is_defined(mt.failing_samples_qc)
        )
        if pheno_col is not None:
            case_filter = sample_filter & (mt[pheno_col] == True) & (hl.is_defined(mt[pheno_col]))
            cont_filter = sample_filter + (mt[pheno_col] == False) & (hl.is_defined(mt[pheno_col]))

            mt = mt.annotate_rows(**{case_het_gt_count: hl.agg.filter(case_filter, hl.agg.count_where(het_gt_filters))})
            mt = mt.annotate_rows(**{cont_het_gt_count: hl.agg.filter(cont_filter, hl.agg.count_where(het_gt_filters))})

        else:
            mt = mt.annotate_rows(**{het_gt_cnt: hl.agg.filter(sample_filter, hl.agg.count_where(het_gt_filters))})
    else:
        mt = mt.annotate_rows(**{het_gt_cnt: hl.agg.count_where(het_gt_filters)})

    ###########################################
    # Define filter conditions for variant ab #
    ###########################################
    het_ab_cond = (
            (((mt.AD[0] / hl.sum(mt.AD)) < min_het_ref_reads) | ((mt.AD[0] / hl.sum(mt.AD)) > max_het_ref_reads))
            & hl.is_defined(mt.AD) & mt.GT.is_het() & hl.is_defined(mt.GT)
    )
    hom_ab_cond = (
            ((mt.AD[0] / hl.sum(mt.AD)) < min_hom_ref_ref_reads) & hl.is_defined(mt.AD) &
            mt.GT.is_hom_ref() & hl.is_defined(mt.GT)
    )
    homalt_ab_cond = (
            ((mt.AD[0] / hl.sum(mt.AD)) > max_hom_alt_ref_reads) & hl.is_defined(mt.AD) &
            mt.GT.is_hom_var() & hl.is_defined(mt.GT)
    )

    #####################################################
    # Get percent of het gts in AB if het gt count != 0 #
    #####################################################
    passing_het_gts = mt.GT.is_het() & hl.is_defined(mt.GT) & (het_ab_cond | hom_ab_cond | homalt_ab_cond)

    if samples_qc:
        sample_filter = ((mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) == 0) &
                         hl.is_defined(mt.pop_outlier_sample) & hl.is_defined(mt.failing_samples_qc))

        if pheno_col is not None:
            case_filter = sample_filter + (mt[pheno_col] == True) & (hl.is_defined(mt[pheno_col]))
            cont_filter = sample_filter + (mt[pheno_col] == False) & (hl.is_defined(mt[pheno_col]))

            mt = mt.annotate_rows(**{
                case_het_gt_ab: hl.agg.filter(case_filter,
                                              hl.cond(
                                                  (mt[case_het_gt_count] > 0) & hl.is_defined(mt[case_het_gt_count]),
                                                  hl.agg.count_where(passing_het_gts) / mt[case_het_gt_count],
                                                  hl.null(hl.tfloat)))})

            mt = mt.annotate_rows(**{
                cont_het_gt_ab: hl.agg.filter(cont_filter,
                                              hl.cond(
                                                  (mt[cont_het_gt_count] > 0) & hl.is_defined(mt[cont_het_gt_count]),
                                                  hl.agg.count_where(passing_het_gts) / mt[cont_het_gt_count],
                                                  hl.null(hl.tfloat)))})

        else:
            mt = mt.annotate_rows(
                **{het_ab: hl.agg.filter(sample_filter,
                                         hl.cond((mt[het_gt_cnt] > 0) & hl.is_defined(mt[het_gt_cnt]),
                                                 hl.agg.count_where(passing_het_gts) / mt[het_gt_cnt],
                                                 hl.null(hl.tfloat)))})
    else:
        mt = mt.annotate_rows(**{het_ab: hl.cond((mt[het_gt_cnt] > 0) & hl.is_defined(mt[het_gt_cnt]),
                                                 hl.agg.count_where(passing_het_gts) / mt[het_gt_cnt],
                                                 hl.null(hl.tfloat))})

    if pheno_col is not None:
        # Check annotations defined for all variants
        undefined_case_het_ct = mt.aggregate_rows(hl.agg.count_where(~(hl.is_defined(mt[case_het_gt_count]))))
        undefined_case_het_ab = mt.aggregate_rows(hl.agg.count_where(~(hl.is_defined(mt[case_het_gt_ab]))))

        if undefined_case_het_ct > 0:
            logging.info(f"Warning- {undefined_case_het_ct} variants have undefined case het GT count.")

        if undefined_case_het_ab > 0:
            logging.info(f"Warning- {undefined_case_het_ab} variants have undefined case het GT allelic balance.")

        undefined_cont_het_ct = mt.aggregate_rows(hl.agg.count_where(~(hl.is_defined(mt[cont_het_gt_count]))))
        undefined_cont_het_ab = mt.aggregate_rows(hl.agg.count_where(~(hl.is_defined(mt[cont_het_gt_ab]))))

        if undefined_cont_het_ct > 0:
            logging.info(f"Warning- {undefined_cont_het_ct} variants have undefined cont het GT count.")

        if undefined_cont_het_ab > 0:
            logging.info(f"Warning- {undefined_cont_het_ab} variants have undefined cont het GT allelic balance.")

    mt = mt.checkpoint(checkpoint_name + "_variant_het_ab_annotated.mt/", overwrite=True)

    return mt


def annotate_variant_het_ab(mt, checkpoint_name, prefix="", samples_qc=False, pheno_col=None, min_het_ref_reads=0.2,
                            max_het_ref_reads=0.8, min_hom_ref_ref_reads=0.9, max_hom_alt_ref_reads=0.1):
    """
    Annotate percentage of het genotypes per variant that are in  allelic balance. That means that the genotype has ref
    reads within the 20-80% range of all reads.

    :param mt: matrix table to annotate
    :param prefix: prefix to add to variant annotations
    :param checkpoint_name: Name of checkpoint to write to, including bucket name or file dir to write to
    :param samples_qc: has samples QC been run, failing samples and population outliers found?
    :param pheno_col: is there a column annotation giving True/False for case status?
    :param min_het_ref_reads: minimum percent reference reads for het genotype
    :param max_het_ref_reads: maximum percent reference reads for het genotype
    :param min_hom_ref_ref_reads: minimum percent reference reads for hom ref genotype
    :param max_hom_alt_ref_reads: maximum percent reference reads for hom alt genotype
    """
    ################################################################
    # Annotate rows with het GT count, then het GT allelic balance #
    ################################################################
    logging.info("Annotating variants with fraction of het genotypes in allelic balance, after filtering out "
                 "population outliers and samples failing samples QC (if run already) and genotypes failing on "
                 "depth and quality by depth measures.")

    mt = count_variant_ab(mt, checkpoint_name, prefix=prefix, samples_qc=samples_qc, pheno_col=None,
                          min_het_ref_reads=min_het_ref_reads, max_het_ref_reads=max_het_ref_reads,
                          min_hom_ref_ref_reads=min_hom_ref_ref_reads, max_hom_alt_ref_reads=max_hom_alt_ref_reads)

    ################################################################################
    # Annotate het GT and het GT allelic balance for cases and controls separately #
    ################################################################################
    if (pheno_col is not None) and (samples_qc is True):
        mt = count_variant_ab(mt, checkpoint_name, prefix=prefix, samples_qc=samples_qc, pheno_col=pheno_col,
                              min_het_ref_reads=min_het_ref_reads, max_het_ref_reads=max_het_ref_reads,
                              min_hom_ref_ref_reads=min_hom_ref_ref_reads, max_hom_alt_ref_reads=max_hom_alt_ref_reads)

    return mt


def find_failing_genotypes_ab(mt, checkpoint_name, prefix="", max_het_ref_reads=0.2, min_het_ref_reads=0.8,
                              min_hom_ref_ref_reads=0.9, max_hom_alt_ref_reads=0.1, count_failing=True):
    """
    Finds genotypes failing allelic balance, defined as the percentage of ref reads that should correspond to that
    called genotype. For hets, percent of ref reads should be between 20-80%, for hom ref they should be at least 90%,
    for hom alt they should be maximum 10%.
    :param mt: matrix table to annotated
    :param args: args to get thresholds from
    :return: 
    """
    if (not prefix.endswith("_")) and (prefix != ""):
        prefix = prefix + "_"
    ###############################
    # Log and annotate thresholds #
    ###############################
    logging.info(f"\nFinding failing het genotypes with ref reads > {max_het_ref_reads}"
                 f"\nor het ref reads < {min_het_ref_reads}"
                 f"\nor failing hom ref genotypes with ref reads < {min_hom_ref_ref_reads}"
                 f"\nor failing hom alt genotypes with ref reads > {max_hom_alt_ref_reads}")

    gt_total = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))

    ############################
    # Define filter conditions #
    ############################
    het_ab_cond = (
            (((mt.AD[0] / hl.sum(mt.AD)) < min_het_ref_reads) | ((mt.AD[0] / hl.sum(mt.AD)) > max_het_ref_reads))
            & hl.is_defined(mt.AD) & mt.GT.is_het() & hl.is_defined(mt.GT)
    )
    hom_ab_cond = (
            ((mt.AD[0] / hl.sum(mt.AD)) < min_hom_ref_ref_reads) & hl.is_defined(mt.AD) &
            mt.GT.is_hom_ref() & hl.is_defined(mt.GT)
    )
    homalt_ab_cond = (
            ((mt.AD[0] / hl.sum(mt.AD)) > max_hom_alt_ref_reads) & hl.is_defined(mt.AD) &
            mt.GT.is_hom_var() & hl.is_defined(mt.GT)
    )

    ###########################
    # Count failing genotypes #
    ###########################
    if count_failing:
        # Get starting genotypes count
        gthet = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_het()))
        gthomvar = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_hom_var()))
        gthomref = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_hom_ref()))

        # Filter het genotypes on allelic balance
        hets_failing_ab = mt.aggregate_entries(hl.agg.count_where(het_ab_cond))
        het_failing_ab_perc = round(hets_failing_ab / gthet * 100, 2)

        # Filter hom ref genotypes on allelic balance
        homref_failing_ab = mt.aggregate_entries(hl.agg.count_where(hom_ab_cond))
        homref_failing_ab_perc = round(homref_failing_ab / gthomref * 100, 4)

        # Filter hom var genotypes on allelic balance
        homalt_failing_ab = mt.aggregate_entries(hl.agg.count_where(homalt_ab_cond))
        homalt_failing_ab_perc = round(homalt_failing_ab / gthomvar * 100, 2)

        # Find genotypes defined but missing AD
        missing_ad_het = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & ~(hl.is_defined(mt.AD)) &
                                                                 mt.GT.is_het()))
        missing_ad_homref = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & ~(hl.is_defined(mt.AD)) &
                                                                    mt.GT.is_hom_ref()))
        missing_ad_homalt = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & ~(hl.is_defined(mt.AD)) &
                                                                    mt.GT.is_hom_var()))

        # Report number and percent of genotypes excluded
        logging.info("\nGenotypes failing ab, after excluding low DP and low GQ genotypes:")
        logging.info(f"Number of het GTs failing ab: {hets_failing_ab}({het_failing_ab_perc}%)")
        logging.info(f"Number of hom ref GTs failing_ab: {homref_failing_ab} ({homref_failing_ab_perc}%)")
        logging.info(f"Number of hom alt GTs failing_ab: {homalt_failing_ab} ({homalt_failing_ab_perc})%")

        logging.info(f"Number of het GTs missing AD info: {missing_ad_het} "
                     f"({round(missing_ad_het/gthet*100, 2)}% of het GTs)")
        logging.info(f"Number of hom ref GTs missing AD info: {missing_ad_homref} "
                     f"({round(missing_ad_homref/gthomref*100, 2)}% of hom ref GTs)")
        logging.info(f"Number of hom alt GTs missing AD info: {missing_ad_homalt} "
                     f"({round(missing_ad_homalt/gthomvar*100, 2)}% of hom alt GTs)")

        global_fail_annot = prefix + "genotype_qc_failing_ab"
        mt = mt.annotate_globals(**{global_fail_annot:
                                   {'het_excluded_ct_percent': [hets_failing_ab, het_failing_ab_perc],
                                    'hom_ref_excluded_ct_percent': [homref_failing_ab, homref_failing_ab_perc],
                                    'hom_var_excluded_ct_percent': [homalt_failing_ab, homalt_failing_ab_perc],
                                    'het_missing_ad': [float(missing_ad_het)],
                                    'hom_ref_missing_ad': [float(missing_ad_homref)],
                                    'hom_alt_missing_ad': [float(missing_ad_homalt)]}})

        #########################################################################
        # re-define filter conditions, global annotation causes source mismatch #
        #########################################################################
        het_ab_cond = (
                (((mt.AD[0] / hl.sum(mt.AD)) < min_het_ref_reads) | ((mt.AD[0] / hl.sum(mt.AD)) > max_het_ref_reads))
                & hl.is_defined(mt.AD) & mt.GT.is_het() & hl.is_defined(mt.GT))
        hom_ab_cond = (
                ((mt.AD[0] / hl.sum(mt.AD)) < min_hom_ref_ref_reads) & hl.is_defined(mt.AD) &
                mt.GT.is_hom_ref() & hl.is_defined(mt.GT))
        homalt_ab_cond = (
                ((mt.AD[0] / hl.sum(mt.AD)) > max_hom_alt_ref_reads) & hl.is_defined(mt.AD) &
                mt.GT.is_hom_var() & hl.is_defined(mt.GT))

    ############################
    # Filter failing genotypes #
    ############################
    mt = mt.filter_entries(het_ab_cond | hom_ab_cond | homalt_ab_cond)
    mt = mt.checkpoint(checkpoint_name + "_GT_ab_filtered.mt/", overwrite=True)

    passing_gts = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    passing_gts_perc = round(passing_gts / gt_total * 100, 2)

    logging.info(f"Number of passing genotypes: {passing_gts} ({passing_gts_perc}%)")

    return mt


def find_failing_vars(mt, checkpoint_name, prefix="", pheno_col=None, count_failing=False, snp_qd=2, indel_qd=3,
                      ab_allowed_dev_het=0.8, p_hwe=1e-6, sex_aware=False, min_call_rate=0.8):
    #######################################################
    # Run hail variant QC method on filtered matrix table #
    #######################################################
    if (not prefix.endswith("_")) and (prefix != ""):
        prefix = prefix + "_"

    failing_name = prefix + "failing_variant_qc"
    varqc_name = prefix + "variant_qc"

    mt = hl.variant_qc(mt, name=varqc_name)

    ######################################################
    # Set HWE calculation in either all or only controls #
    ######################################################
    if pheno_col is not None:
        case_filters = (mt[pheno_col] == False) & hl.is_defined(mt[pheno_col])
        mt = mt.annotate_rows(
            hwe_ctrls_only=hl.agg.filter(case_filters, hl.agg.hardy_weinberg_test(mt.GT)))

        num_controls = mt.aggregate_cols(hl.agg.count_where((mt[pheno_col] == False) & hl.is_defined(mt[pheno_col])))

        logging.info(f'Calculated HWE with {num_controls} controls only '
                     f'(ignoring samples with missing phenotype information)')
    else:
        logging.info('Phenotype column name not provided, calculating HWE for variants using all samples.')

    # Instantiate annotation, count starting variants
    mt = mt.annotate_rows(**{failing_name: hl.empty_array(hl.tstr)})
    total_variants = mt.count_rows()

    ########################
    # Failing VQSR filters #
    ########################
    vqsr_defined = mt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt.filters)))
    if vqsr_defined == 0:
        logging.error(
            f"Note! mt.filters annotation undefined for all variants! Not filtering variants on this measure.")
    else:
        vqsr_filter = hl.cond((hl.len(mt.filters) != 0) & hl.is_defined(mt.filters),
                              mt[failing_name].append("failing_VQSR_filters"),
                              mt[failing_name])
        mt = mt.annotate_rows(**{failing_name: vqsr_filter})

        vqsr_miss = hl.cond(~hl.is_defined(mt.filters), mt[failing_name].append("missing_VQSR_filters"),
                            mt[failing_name])
        mt = mt.annotate_rows(**{failing_name: vqsr_miss})

    #####################################################################
    # Find variants that are snps or indels with QD less than threshold #
    #####################################################################
    snps = hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.info.QD < snp_qd) & hl.is_defined(mt.info.QD)
    indels = (~hl.is_snp(mt.alleles[0], mt.alleles[1])) & (mt.info.QD < indel_qd) & hl.is_defined(mt.info.QD)

    qd_filter = hl.cond(snps | indels, mt[failing_name].append("failing_QD"), mt[failing_name])
    mt = mt.annotate_rows(**{failing_name: qd_filter})

    ###################################################################
    # Find variants with >20% of het genotypes out of allelic balance #
    ###################################################################
    failing_het_ab_name = prefix + 'frac_het_gts_in_ab'
    ab_filter = hl.cond(hl.is_defined(mt[failing_het_ab_name]) &
                        (mt[failing_het_ab_name] < ab_allowed_dev_het),
                        mt[failing_name].append("failing_het_ab"), mt[failing_name])
    mt = mt.annotate_rows(**{failing_name: ab_filter})

    ####################
    # Failing HWE test #
    ####################
    if pheno_col is not None:
        hwe_cond = (mt.hwe_ctrls_only.p_value < p_hwe) & hl.is_defined(mt.hwe_ctrls_only.p_value)
    else:
        hwe_cond = (mt.row[varqc_name].p_value_hwe < p_hwe) & hl.is_defined(mt[varqc_name].p_value_hwe)

    hwe_filter = hl.cond(hwe_cond, mt[failing_name].append("failing_hwe"), mt[failing_name])
    mt = mt.annotate_rows(**{failing_name: hwe_filter})

    ##############################################
    # Find variants not passing call rate filter #
    ##############################################
    if sex_aware:
        call_rate_cond = (mt.sexaware_call_rate < min_call_rate) & hl.is_defined(mt.sexaware_call_rate)
    else:
        call_rate_cond = (mt[varqc_name].call_rate < min_call_rate) & hl.is_defined(mt[varqc_name].call_rate)

    call_rate_filter = hl.cond(call_rate_cond, mt[failing_name].append("failing_call_rate"), mt[failing_name])
    mt = mt.annotate_rows(**{failing_name: call_rate_filter})

    ##############
    # Checkpoint #
    ##############
    mt = mt.checkpoint(checkpoint_name + "_failing_vars_annotated.mt/", overwrite=True)

    ##########################
    # Count failing variants #
    ##########################
    if count_failing:
        filter_dict = {}

        # VQSR Filters
        failing_vqsr = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains("failing_VQSR_filters")))
        filter_dict["failing_VQSR"] = failing_vqsr
        failing_vqsr_perc = round(failing_vqsr / total_variants * 100, 2)

        missing_vqsr = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains("missing_VQSR_filters")))
        filter_dict['missing_VQSR'] = missing_vqsr
        missing_vqsr_perc = round(missing_vqsr / total_variants * 100, 2)

        # QD Filters
        failing_qd = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains("failing_QD")))
        failing_qd_perc = round(failing_qd / total_variants * 100, 2)
        filter_dict["failing_QD"] = failing_qd

        qd_defined = mt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt.info.QD)))
        if qd_defined < total_variants:
            qd_miss = hl.cond(~(hl.is_defined(mt.info.QD)), mt[failing_name].append("missing_QD"),
                              mt[failing_name])

            mt = mt.annotate_rows(**{failing_name: qd_miss})

            missing_qd = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains("missing_QD")))
            filter_dict['missing_QD'] = missing_qd
            missing_qd_perc = round(missing_qd / total_variants * 100, 2)

            logging.info(f"Note! QD missing for some variants: {missing_qd}, {missing_qd_perc}%")

        # allelic balance filters
        failing_ab = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains("failing_het_ab")))
        failing_ab_perc = round(failing_ab / total_variants * 100, 2)
        filter_dict["failing_het_allelic_balance"] = failing_ab

        ab_defined = mt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt[failing_het_ab_name])))
        if ab_defined < total_variants:
            logging.info(f"Note! mt.{failing_het_ab_name} annotation defined for only {ab_defined} variants! "
                          f"Variants missing this annotation not filtered on this measure.")

        # Hardy-Weinberg equilibrium filters
        if pheno_col is not None:
            hwe_defined = mt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt.hwe_ctrls_only.p_value)))
        else:
            hwe_defined = mt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt[varqc_name].p_value_hwe)))

        failing_hwe = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains('failing_hwe')))
        failing_hwe_perc = round(failing_hwe / total_variants * 100, 2)
        filter_dict["failing_hwe"] = failing_hwe

        if hwe_defined < total_variants:
            logging.error(f"Note! HWE annotation defined for only {hwe_defined} variants! "
                          f"Something is wrong, check this.")

        # Call rate
        if sex_aware:
            cr_defined = mt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt.sexaware_call_rate)))
            low_passing_vars = mt.aggregate_rows(hl.agg.count_where(hl.len(mt.low_pass_failing_variant_qc) == 0))
            if cr_defined < low_passing_vars:
                logging.error(f"Note! sex aware call rate annotation defined for only {cr_defined} variants, where "
                              f"there are {low_passing_vars} variants passing low pass variant QC. Something is wrong!")
        else:
            cr_defined = mt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt[varqc_name].call_rate)))
            if cr_defined < total_variants:
                logging.error(f"Note! Call rate annotation defined for only {cr_defined} variants! "
                              f"Something is wrong, check this.")

        failing_call_rate = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains("failing_call_rate")))
        failing_cr_perc = round(failing_call_rate / total_variants * 100, 2)
        filter_dict["failing_call_rate"] = failing_call_rate

        if pheno_col is not None:
            hwe_tag = "controls samples only"
        else:
            hwe_tag = "all samples"

        mt = mt.annotate_globals(**{prefix + "failing_variant_qc_counts": filter_dict})

        logging.info(f"\nVariants failing VQSR filters: {failing_vqsr}, {failing_vqsr_perc}%"
                     f"\nsnp QD < {snp_qd} or indel QD < {indel_qd}: {failing_qd}, {failing_qd_perc}%"
                     f"\n>{ab_allowed_dev_het * 100}% of het genotypes out of allelic balance: "
                     f"{failing_ab}, {failing_ab_perc}%"
                     f"\np value HWE < {p_hwe} in {hwe_tag}: {failing_hwe}, {failing_hwe_perc}%"
                     f"\ncall rate < {min_call_rate}: {failing_call_rate}, {failing_cr_perc}%"
                     f"\n\nVariants missing VQSR filters: {missing_vqsr}, {missing_vqsr_perc}%"
                     f"\n")

        passing = mt.aggregate_rows(hl.agg.count_where(hl.len(mt[failing_name]) == 0))

        logging.info(f"Variants passing QC: {passing}, {round(passing / total_variants * 100, 2)}%")

    return mt


def variant_quality_control(
        mt, checkpoint_name, annotation_prefix="", min_dp=10, min_gq=20,
        max_het_ref_reads=0.8, min_het_ref_reads=0.2, min_hom_ref_ref_reads=0.9, max_hom_alt_ref_reads=0.1,
        call_rate=0.8, p_hwe=1e-6, snp_qd=2, indel_qd=3, ab_allowed_dev_het=0.8,filter_missing_measures=False,
        count_failing=True, sex_aware_call_rate=False, pheno_col=None, samples_qc=False, force=False):
    """
    Function to find variants failing on QC measures, which can be run in 'low_pass' or 'final' mode, with varying
    filters given by args depending on the mode.

    :param mt: matrix table to filter
    :param args: arguments object with thresholds
    :return: returns matrix table with entries and variants annotated with info on whether they fail QC metrics.
    """
    ############################################################################
    # Check that multiallelic variants have been split, define annotation name #
    ############################################################################
    utils.check_multi_split(mt)

    if annotation_prefix is not None:
        failing_name = annotation_prefix + "_failing_variant_qc"
    else:
        failing_name = "failing_variant_qc"

    if (not annotation_prefix.endswith("_")) and (annotation_prefix != ""):
        annotation_prefix = annotation_prefix + "_"

    if checkpoint_name.endswith(".mt"):
        checkpoint_name = checkpoint_name.replace(".mt", "")
    elif checkpoint_name.endswith(".mt/"):
        checkpoint_name = checkpoint_name.replace(".mt/", "")

    #####################################
    # Annotate QC thresholds to globals #
    #####################################
    if pheno_col is not None:
        hwe_tag = "controls samples only"
    else:
        hwe_tag = "all samples"

    var_thresholds = {'pval_hwe': str(p_hwe), 'hwe_excluded_in': str(hwe_tag),
                      'call_rate': str(call_rate),
                      'sex_aware_call_rate': str(sex_aware_call_rate), 'snp_qd': str(snp_qd),
                      'indel_qd': str(indel_qd),
                      'het_max_ref_reads_thresh': str(max_het_ref_reads),
                      'het_min_ref_reads_thresh': str(min_het_ref_reads),
                      'perc_het_reads_ab_allowed': str(ab_allowed_dev_het)}

    mt = mt.annotate_globals(**{annotation_prefix + "variant_qc_thresholds": hl.literal(var_thresholds)})

    gt_thresholds = {'min_DP': min_dp, 'min_GQ': min_gq, 'max_het_ref_reads': max_het_ref_reads,
                     'min_het_ref_reads': min_het_ref_reads, 'min_hom_ref_ref_reads': min_hom_ref_ref_reads,
                     'max_hom_alt_ref_reads': max_hom_alt_ref_reads }

    mt = mt.annotate_globals(**{annotation_prefix + "genotype_qc_thresholds": hl.literal(gt_thresholds)})

    ############################################################
    # Find failing genotypes and do allelic balance annotation #
    ############################################################
    # Filter genotypes failing on depth + quality
    if (not utils.check_exists(checkpoint_name + "_DP_GQ_filtered.mt/")) or force:
        mt_gtfilt = filter_failing_GTs_depth_quality(
            mt, checkpoint_name, prefix=annotation_prefix, min_dp=min_dp, min_gq=min_gq,
            filter_missing_measures=filter_missing_measures, count_failing=count_failing)
    else:
        logging.info("Detected DP + GQ filtered mt exists, loading that.")
        mt_gtfilt = hl.read_matrix_table(checkpoint_name + "_DP_GQ_filtered.mt/")

    # Annotate variants failing het AB measure (percentage of het GT calls for that variant *in balance*)
    if (not utils.check_exists(checkpoint_name + "_variant_het_ab_annotated.mt/")) or force:
        mt_abannot = annotate_variant_het_ab(
            mt_gtfilt, checkpoint_name, prefix=annotation_prefix, samples_qc=samples_qc, pheno_col=pheno_col,
            max_het_ref_reads=max_het_ref_reads, min_het_ref_reads=min_het_ref_reads,
            min_hom_ref_ref_reads=min_hom_ref_ref_reads, max_hom_alt_ref_reads=max_hom_alt_ref_reads
        )
    else:
        logging.info("Detected het AB annotated mt exists, loading that.")
        mt_abannot = hl.read_matrix_table(checkpoint_name + "_variant_het_ab_annotated.mt/")

    # Filter genotypes failing on allelic balance
    if (not utils.check_exists(checkpoint_name + "_GT_ab_filtered.mt/")) or force:
        mt_gtfilt2 = find_failing_genotypes_ab(
            mt_abannot, checkpoint_name, prefix=annotation_prefix, max_het_ref_reads=max_het_ref_reads,
            min_het_ref_reads=min_het_ref_reads, min_hom_ref_ref_reads=min_hom_ref_ref_reads,
            max_hom_alt_ref_reads=max_hom_alt_ref_reads, count_failing=count_failing
        )
    else:
        logging.info("Detected genotype AB filtered mt exists, loading that.")
        mt_gtfilt2 = hl.read_matrix_table(checkpoint_name + "_GT_ab_filtered.mt/")

    ###############################################
    # Filter matrix table to only passing samples #
    ###############################################
    if samples_qc:
        logging.info("Filtering out samples failing QC.")
        mt_gtfilt2 = sq.filter_failing(
            mt_gtfilt2, checkpoint_name, prefix=annotation_prefix, pheno_col=pheno_col, entries=False, variants=False,
            samples=True, pheno_qc=False, min_dp=min_dp, max_het_ref_reads=max_het_ref_reads,
            min_het_ref_reads=min_het_ref_reads, min_hom_ref_ref_reads=min_hom_ref_ref_reads,
            max_hom_alt_ref_reads=max_hom_alt_ref_reads)

    #####################################
    # Annotate rest of failing variants #
    #####################################
    if (not utils.check_exists(checkpoint_name + "_failing_vars_annotated.mt/")) or force:
        mt_varannot = find_failing_vars(mt_gtfilt2, checkpoint_name, prefix=annotation_prefix, pheno_col=pheno_col,
                                        count_failing=count_failing, snp_qd=snp_qd, indel_qd=indel_qd,
                                        ab_allowed_dev_het=ab_allowed_dev_het, p_hwe=p_hwe, sex_aware=sex_aware_call_rate,
                                        min_call_rate=call_rate)
    else:
        logging.info("Detected failing variant mt exists, loading that.")
        mt_varannot = hl.read_matrix_table(checkpoint_name + "_failing_vars_annotated.mt/")

    #########################################################################
    # Add annotations back to main matrix table (not filtered), return that #
    #########################################################################
    varqc_name = annotation_prefix + "variant_qc"
    mt = mt.annotate_rows(**{varqc_name: mt_varannot.index_rows(mt.row_key)[varqc_name]})
    mt = mt.annotate_rows(**{failing_name: mt_varannot.index_rows(mt.row_key)[failing_name]})

    if (pheno_col is not None) and (samples_qc is True):
        case_het_gt_ab = annotation_prefix + 'case_frac_het_gts_in_ab'
        mt = mt.annotate_rows(**{case_het_gt_ab: mt_varannot.index_rows(mt.row_key)[case_het_gt_ab]})
        cont_het_gt_ab = annotation_prefix + 'control_frac_het_gts_in_ab'
        mt = mt.annotate_rows(**{cont_het_gt_ab: mt_varannot.index_rows(mt.row_key)[cont_het_gt_ab]})

    else:
        het_ab = annotation_prefix + 'frac_het_gts_in_ab'
        mt = mt.annotate_rows(**{het_ab: mt_varannot.index_rows(mt.row_key)[het_ab]})

    if count_failing:
        dp_qc_name = annotation_prefix + "genotype_qc_failing_quality_depth"
        mt = mt.annotate_globals(**{dp_qc_name: mt_varannot.index_globals()[dp_qc_name]})
        ab_name = annotation_prefix + "genotype_qc_failing_ab"
        mt = mt.annotate_globals(**{ab_name: mt_varannot.index_globals()[ab_name]})
        var_name = annotation_prefix + "failing_variant_qc_counts"
        mt = mt.annotate_globals(**{var_name: mt_varannot.index_globals()[var_name]})

    if pheno_col is not None:
        mt = mt.annotate_rows(hwe_ctrls_only=mt_varannot.index_rows(mt.row_key).hwe_ctrls_only)

    return mt


def find_variants_failing_by_pheno(mt, args):
    """
    Finds variants failing call rate filters specifically by phenotype.
    :param mt:
    :param args:
    :return:
    """
    # Check that samples QC has been run
    try:
        test = hl.is_defined(mt.failing_samples_qc)
        test = hl.is_defined(mt.pop_outlier_sample)
    except Exception as e:
        logging.info('failing_samples_qc and pop_outlier_sample not defined! Run samples QC before running this.')
        logging.info(e)

    # Initialize sample annotation for phenotype-specific measures
    mt = mt.annotate_rows(failing_pheno_varqc=hl.empty_array(hl.tstr))
    # Count total variants
    total_vars = mt.count_rows()

    ##############################################################
    # Annotate variants failing on case-specific allelic balance #
    ##############################################################
    case_cond = hl.cond(hl.is_defined(mt.final_case_frac_het_gts_in_ab) &
                        (mt.final_case_frac_het_gts_in_ab < args.ab_allowed_dev_het),
                        mt.failing_pheno_varqc.append("failing_case_het_ab"),
                        mt.failing_pheno_varqc)
    mt = mt.annotate_rows(failing_pheno_varqc=case_cond)
    failing_ab_case = mt.aggregate_rows(hl.agg.count_where(mt.failing_pheno_varqc.contains("failing_case_het_ab")))
    ab_case_perc = round(failing_ab_case/total_vars*100, 2)

    control_cond = hl.cond(hl.is_defined(mt.final_control_frac_het_gts_in_ab) &
                           (mt.final_control_frac_het_gts_in_ab < args.ab_allowed_dev_het),
                           mt.failing_pheno_varqc.append("failing_control_het_ab"),
                           mt.failing_pheno_varqc)
    mt = mt.annotate_rows(failing_pheno_varqc=control_cond)
    failing_ab_cont = mt.aggregate_rows(hl.agg.count_where(mt.failing_pheno_varqc.contains("failing_control_het_ab")))
    ab_case_cont = round(failing_ab_cont/total_vars*100, 2)

    ################################################################
    # Find variants failing on case-specific call rate (sex aware) #
    ################################################################
    case_cr_cond = hl.cond(hl.is_defined(mt.sexaware_case_call_rate) &
                           (mt.sexaware_case_call_rate < args.pheno_call_rate),
                           mt.failing_pheno_varqc.append("failing_case_call_rate"),
                           mt.failing_pheno_varqc)
    mt = mt.annotate_rows(failing_pheno_varqc=case_cr_cond)
    failing_cr_case = mt.aggregate_rows(hl.agg.count_where(mt.failing_pheno_varqc.contains("failing_case_call_rate")))
    cr_case_perc = round(failing_cr_case/total_vars*100, 2)

    cont_cr_cond = hl.cond(hl.is_defined(mt.sexaware_cont_call_rate) &
                           (mt.sexaware_cont_call_rate < args.pheno_call_rate),
                           mt.failing_pheno_varqc.append("failing_cont_call_rate"),
                           mt.failing_pheno_varqc)
    mt = mt.annotate_rows(failing_pheno_varqc=cont_cr_cond)
    failing_cr_cont = mt.aggregate_rows(hl.agg.count_where(mt.failing_pheno_varqc.contains("failing_cont_call_rate")))
    cr_cont_perc = round(failing_cr_cont/total_vars*100, 2)

    mt = mt.annotate_globals(case_control_callrate_threshold=args.pheno_call_rate)

    logging.info(f"Number of variants failing on case ab: {failing_ab_case} ({ab_case_perc}%)")
    logging.info(f"Number of variants failing on control ab: {failing_ab_cont} ({ab_case_cont}%)")
    logging.info(f"Number of variants failing on case call rate: {failing_cr_case} ({cr_case_perc}%)")
    logging.info(f"Number of variants failing on control call rate: {failing_cr_cont} ({cr_cont_perc}%)")

    failing_any = mt.aggregate_rows(hl.agg.count_where(hl.len(mt.failing_pheno_varqc) != 0))
    failing_perc = round(failing_any/total_vars*100, 2)
    logging.info(f"Number of variants failing on any phenotype-specific measure: {failing_any} ({failing_perc}%)")

    return mt
