"""
Functions for variant QC for exome sequencing data with Hail.

Author: Lea Urpa, August 2020
"""
import logging
import hail as hl
import samples_qc as sq


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

    mt_ldpruned = mt_ldpruned.annotate_globals(ld_pruning_parameters={'r2': args.r2,
                                                                      'bp_window_size': args.bp_window_size})
    return mt_ldpruned


def find_failing_genotypes_depth_quality(mt, args, prefix):
    """
    Annotates genotypes for initial low-pass variant QC in pipeline. Finds genotypes with low depth and GQ values.
    :param mt: matrix table to filter
    :return: returns filtered matrix table
    """
    # Log thresholds + report
    logging.info(f"Finding genotypes with min_dp < {args.min_dp}, or min GQ < {args.min_gq}")
    globals_annot_name = prefix + "_genotype_qc_thresh_dp_gq"
    failing_name = prefix + "_failing_depth_quality"
    mt = mt.annotate_globals(**{globals_annot_name: {'min_dp': args.min_dp, 'min_gq': args.min_gq}})

    # Get starting genotypes count, instantiate genotype annotation
    gt_total = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    gt_het_homalt = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) &
                                                            (mt.GT.is_hom_var() | mt.GT.is_het())))
    gt_homref = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_hom_ref()))
    mt = mt.annotate_entries(**{failing_name: hl.empty_array(hl.tstr)})

    ############################
    # Filter entries for depth #
    ############################
    depth_cond = hl.cond((mt.DP < args.min_dp) & hl.is_defined(mt.DP),
                         mt[failing_name].append("failing_DP"), mt[failing_name])
    mt = mt.annotate_entries(**{failing_name: depth_cond})
    missing_dp = hl.cond(~(hl.is_defined(mt.DP)) & hl.is_defined(mt.GT),
                         mt[failing_name].append("missing_DP"), mt[failing_name])
    mt = mt.annotate_entries(**{failing_name: missing_dp})

    failing_depth = mt.aggregate_entries(hl.agg.count_where(mt[failing_name].contains("failing_DP")))
    failing_depth_perc = round(failing_depth / gt_total*100, 2)
    missing_depth = mt.aggregate_entries(hl.agg.count_where(mt[failing_name].contains("missing_DP")))
    missing_depth_perc = round(missing_depth / gt_total*100, 2)

    ###########################################
    # Filter GQ or PL for different genotypes #
    ###########################################
    gq_cond = hl.cond((mt.GT.is_het() | mt.GT.is_hom_var()) & (mt.PL[0] < args.min_gq) &
                      hl.is_defined(mt.GT) & hl.is_defined(mt.PL),
                      mt[failing_name].append("failing_PL"), mt[failing_name])
    mt = mt.annotate_entries(**{failing_name: gq_cond})

    pl_cond = hl.cond(mt.GT.is_hom_ref() & (mt.GQ < args.min_gq) & hl.is_defined(mt.GQ) & hl.is_defined(mt.GT),
                      mt[failing_name].append("failing_GQ"), mt[failing_name])
    mt = mt.annotate_entries(**{failing_name: pl_cond})

    gq_miss = hl.cond(~(hl.is_defined(mt.GQ)) & hl.is_defined(mt.GT) & mt.GT.is_hom_ref(),
                      mt[failing_name].append("missing_GQ_hom_ref"), mt[failing_name])

    mt = mt.annotate_entries(**{failing_name: gq_miss})

    pl_miss = hl.cond(~(hl.is_defined(mt.PL)) & (mt.GT.is_het() | mt.GT.is_hom_var()) & hl.is_defined(mt.GT),
                      mt[failing_name].append("missing_PL_het_homvar"), mt[failing_name])
    mt = mt.annotate_entries(**{failing_name: pl_miss})

    failing_gq = mt.aggregate_entries(hl.agg.count_where(mt[failing_name].contains("failing_GQ")))
    failing_pl = mt.aggregate_entries(hl.agg.count_where(mt[failing_name].contains("failing_PL")))
    missing_gq = mt.aggregate_entries(hl.agg.count_where(mt[failing_name].contains("missing_GQ_or_PL")))
    missing_pl = mt.aggregate_entries(hl.agg.count_where(mt[failing_name].contains("missing_PL_het_homvar")))

    logging.info(f"Number of GTs with depth < {args.min_dp}: {failing_depth} ({failing_depth_perc}%)")
    logging.info(f"Number of hom ref GTs with GQ < {args.min_gq}: {failing_gq} "
                 f"({round(failing_gq/gt_total*100, 2)}% of all GTs, "
                 f"{round(failing_gq/gt_homref*100, 2)}% of homref GTs)")
    logging.info(f"Number of het or hom alt GTs with PL[0] < {args.min_gq}: {failing_pl} "
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

    passing_gts = mt.aggregate_entries(hl.agg.count_where((hl.len(mt[failing_name]) == 0) & hl.is_defined(mt.GT)))
    passing_gts_perc = round(passing_gts / gt_total * 100, 2)

    logging.info(f"Number of passing genotypes: {passing_gts} ({passing_gts_perc}%)")

    globals_fail_annot = prefix + "_genotype_qc_failing_quality_depth"
    mt = mt.annotate_globals(**{globals_fail_annot: {'depth_excluded': failing_depth, 'quality_excluded': failing_gq}})

    return mt


def find_failing_genotypes_ab(mt, args, prefix):
    """
    Finds genotypes failing allelic balance, defined as the percentage of ref reads that should correspond to that
    called genotype. For hets, percent of ref reads should be between 20-80%, for hom ref they should be at least 90%,
    for hom alt they should be maximum 10%.
    :param mt: matrix table to annotated
    :param args: args to get thresholds from
    :return: 
    """
    ###############################
    # Log and annotate thresholds #
    ###############################
    logging.info(f"\nFinding failing het genotypes with ref reads > {args.max_het_ref_reads}"
                 f"\nor het ref reads < {args.min_het_ref_reads}"
                 f"\nor failing hom ref genotypes with ref reads < {args.min_hom_ref_ref_reads}"
                 f"\nor failing hom alt genotypes with ref reads > {args.max_hom_alt_ref_reads}")

    globals_annot_name = prefix + "_genotype_qc_thresholds_ab"
    mt = mt.annotate_globals(**{globals_annot_name: {'max_het_ref_reads': args.max_het_ref_reads,
                                                     'min_het_ref_reads': args.min_het_ref_reads,
                                                     'min_hom_ref_ref_reads': args.min_hom_ref_ref_reads,
                                                     'max_hom_alt_ref_reads': args.max_hom_alt_ref_reads}})

    ##################################################################
    # Check that depth and gq filters have been applied, then filter #
    ##################################################################
    failing_dp_name = prefix + "_failing_depth_quality"
    try:
        test = hl.is_defined(mt[failing_dp_name])
    except Exception as e:
        logging.error("Error! Run find_failing_genotypes_depth_quality() on this matrix table before running "
                      "find_failing_genotypes_ab().")
        logging.error(e)
        exit()

    ################################
    # Get starting genotypes count #
    ################################
    failing_name = prefix + "_failing_ab"
    mt = mt.annotate_entries(**{failing_name: hl.empty_array(hl.tstr)})
    gthet = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_het()))
    gthomvar = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_hom_var()))
    gthomref = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_hom_ref()))

    ###########################################
    # Filter het genotypes on allelic balance #
    ###########################################
    het_ab_cond = (((mt.AD[0] / hl.sum(mt.AD)) < args.min_het_ref_reads) |
                  ((mt.AD[0] / hl.sum(mt.AD)) > args.max_het_ref_reads)) & hl.is_defined(mt.AD)
    ab_cond_het = hl.cond(mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt[failing_dp_name]) == 0) & het_ab_cond,
                          mt[failing_name].append("failing_het_ab"), mt[failing_name])

    mt = mt.annotate_entries(**{failing_name: ab_cond_het})

    hets_failing_ab = mt.aggregate_entries(hl.agg.count_where(mt[failing_name].contains("failing_het_ab")))
    het_failing_ab_perc = round(hets_failing_ab / gthet * 100, 2)

    ###############################################
    # Filter hom ref genotypes on allelic balance #
    ###############################################
    hom_ab_cond = ((mt.AD[0] / hl.sum(mt.AD)) < args.min_hom_ref_ref_reads) & hl.is_defined(mt.AD)
    ab_cond_homref = hl.cond(mt.GT.is_hom_ref() & hom_ab_cond & hl.is_defined(mt.GT),
                             mt[failing_name].append("failing_homref_ab"), mt[failing_name])

    mt = mt.annotate_entries(**{failing_name: ab_cond_homref})

    homref_failing_ab = mt.aggregate_entries(hl.agg.count_where(mt[failing_name].contains("failing_homref_ab")))
    homref_failing_ab_perc = round(homref_failing_ab / gthomref * 100, 2)

    ###############################################
    # Filter hom var genotypes on allelic balance #
    ###############################################
    homalt_ab_cond = ((mt.AD[0] / hl.sum(mt.AD)) > args.max_hom_alt_ref_reads) & hl.is_defined(mt.AD)
    ab_cond_homalt = hl.cond(mt.GT.is_hom_var() & homalt_ab_cond & hl.is_defined(mt.GT),
                             mt[failing_name].append("failing_homalt_ab"), mt[failing_name])

    mt = mt.annotate_entries(**{failing_name: ab_cond_homalt})

    homalt_failing_ab = mt.aggregate_entries(hl.agg.count_where(mt[failing_name].contains("failing_homalt_ab")))
    homalt_failing_ab_perc = round(homalt_failing_ab / gthomvar * 100, 2)

    #########################################
    # Find genotypes defined but missing AD #
    #########################################
    missing_ad_het = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & ~(hl.is_defined(mt.AD)) &
                                                             mt.GT.is_het()))
    missing_ad_homref = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & ~(hl.is_defined(mt.AD)) &
                                                                mt.GT.is_hom_ref()))
    missing_ad_homalt = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & ~(hl.is_defined(mt.AD)) &
                                                                mt.GT.is_hom_var()))

    ###################################################
    # Report number and percent of genotypes excluded #
    ###################################################
    logging.info("\nGenotypes failing ab, after excluding low DP and low GQ genotypes:")
    logging.info(f"Number of het GTs failing ab: {hets_failing_ab}({het_failing_ab_perc}%)")
    logging.info(f"Number of hom ref GTs failing_ab: {homref_failing_ab} ({homref_failing_ab_perc}%)")
    logging.info(f"Number of hom alt GTs failing_ab: {homalt_failing_ab} ({homalt_failing_ab_perc})")

    logging.info(f"Number of het GTs missing AD info: {missing_ad_het} "
                 f"({round(missing_ad_het/gthet*100, 2)}% of het GTs)")
    logging.info(f"Number of hom ref GTs missing AD info: {missing_ad_homref} "
                 f"({round(missing_ad_homref/gthomref*100, 2)}% of hom ref GTs)")
    logging.info(f"Number of hom alt GTs missing AD info: {missing_ad_homalt} "
                 f"({round(missing_ad_homalt/gthomvar*100, 2)}% of hom alt GTs)")

    global_fail_annot = prefix + "_genotype_qc_failing_ab"
    mt = mt.annotate_globals(**{global_fail_annot:
                               {'het_excluded_ct_percent': [hets_failing_ab, het_failing_ab_perc],
                                'hom_ref_excluded_ct_percent': [homref_failing_ab, homref_failing_ab_perc],
                                'hom_var_excluded_ct_percent': [homalt_failing_ab, homalt_failing_ab_perc],
                                'het_missing_ad': [missing_ad_het],
                                'hom_ref_missing_ad': [missing_ad_homref],
                                'hom_var_missing_ad': [missing_ad_homalt]}})

    return mt


def annotate_variant_het_ab(mt, args, prefix):
    """
     Annotate percentage of het genotypes per variant that are in  allelic balance. That means that the genotype has ref
    reads within the 20-80% range of all reads.

    :param mt: matrix table to annotate
    :param prefix: prefix to add to variant annotation
    :return: returns annotated matrix table
    """
    ##########################################
    # Check dataset for expected annotations #
    ##########################################
    logging.info("Annotating variants with fraction of het genotypes in allelic balance, after filtering out "
                 "population outliers and samples failing samples QC (if run already) and genotypes failing on "
                 "depth and quality by depth measures.")
    het_gt_cnt = prefix + '_het_gt_count'
    het_ab = prefix + '_frac_het_gts_in_ab'
    failing_dp_name = prefix + "_failing_depth_quality"
    failing_ab_name = prefix + "_failing_ab"

    # Check that genotype filter annotations have been added
    try:
        test = hl.is_defined(mt[failing_ab_name])
        test = hl.is_defined(mt[failing_dp_name])
    except Exception as e:
        logging.error("Error! find_failing_genotypes_ab() and find_failing_genotypes_depth_quality() should be run "
                      "before annotating variant's % of hets in allelic balance, to filter out bad genotypes.")
        logging.error(e)
        exit()

    # Check if samples QC has been run yet, if not add empty annotations
    try:
        test = hl.is_defined(mt.pop_outlier_sample)
    except Exception as e:
        logging.info("Detected samples QC not run, adding False pop_outlier_sample annotation to columns.")
        logging.info(e)
        mt = mt.annotate_cols(pop_outlier_sample=False)

    try:
        test = hl.is_defined(mt.failing_samples_qc)
    except Exception as e:
        logging.info("Detected samples QC not run, adding empty array failing_samples_qc annotation to columns.")
        logging.info(e)
        mt = mt.annotate_cols(failing_samples_qc=hl.empty_array(hl.tstr))

    ################################################################
    # Annotate rows with het GT count, then het GT allelic balance #
    ################################################################
    # Get count of het genotypes per variant
    case_filter = (mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) == 0) & \
                  hl.is_defined(mt.pop_outlier_sample) & hl.is_defined(mt.failing_samples_qc)
    het_gt_filters = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt[failing_dp_name]) == 0) & \
                     hl.is_defined(mt[failing_dp_name])

    mt = mt.annotate_rows(**{het_gt_cnt: hl.agg.filter(case_filter, hl.agg.count_where(het_gt_filters))})

    # Get percent of het gts in AB if het gt count != 0
    case_filter = (mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) == 0) & \
                  hl.is_defined(mt.pop_outlier_sample) & hl.is_defined(mt.failing_samples_qc)
    passing_het_gts = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt[failing_dp_name]) == 0) & \
                      ~(mt[failing_ab_name].contains("failing_het_ab")) & hl.is_defined(mt[failing_dp_name])

    # Counting het GTs skips variants where there are no het GTs (which leads to div by 0 error)
    mt = mt.annotate_rows(**{het_ab:
             hl.agg.filter(case_filter,
                           hl.cond((mt[het_gt_cnt] > 0) & hl.is_defined(mt[het_gt_cnt]),
                                   hl.agg.count_where(passing_het_gts) / mt[het_gt_cnt],
                                   hl.null(hl.tfloat)))})

    # Check annotations defined for all variants
    undefined_het_ct = mt.aggregate_rows(hl.agg.count_where(~(hl.is_defined(mt[het_gt_cnt]))))
    undefined_het_ab = mt.aggregate_rows(hl.agg.count_where(~(hl.is_defined(mt[het_ab]))))

    if undefined_het_ct > 0:
        logging.info(f"Warning- {undefined_het_ct} variants have undefined het GT count.")

    if undefined_het_ab > 0:
        logging.info(f"Warning- {undefined_het_ab} variants have undefined het GT allelic balance.")

    ################################################################################
    # Annotate het GT and het GT allelic balance for cases and controls separately #
    ################################################################################
    if args.pheno_col is not None:
        case_het_gt_count = prefix + '_case_het_gt_count'
        case_het_gt_ab = prefix + '_case_frac_het_gts_in_ab'
        cont_het_gt_count = prefix + '_control_het_gt_count'
        cont_het_gt_ab = prefix + '_control_frac_het_gts_in_ab'

        # Count case het GTs
        case_filter = (mt[args.pheno_col] == True) & (hl.is_defined(mt[args.pheno_col])) & \
                      (mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) == 0) & \
                      hl.is_defined(mt.pop_outlier_sample) & hl.is_defined(mt.failing_samples_qc)
        het_gt_filters = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt[failing_dp_name]) == 0) & \
                         hl.is_defined(mt[failing_dp_name])

        mt = mt.annotate_rows(**{case_het_gt_count: hl.agg.filter(case_filter, hl.agg.count_where(het_gt_filters))})

        # Get case het GT fraction in allelic balance
        case_filter = (mt[args.pheno_col] == True) & (hl.is_defined(mt[args.pheno_col])) & \
                      (mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) > 0) & \
                      hl.is_defined(mt.pop_outlier_sample) & hl.is_defined(mt.failing_samples_qc)
        passing_het_gts = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt[failing_dp_name]) == 0) & \
                          ~(mt[failing_ab_name].contains("failing_het_ab")) & hl.is_defined(mt[failing_dp_name])

        mt = mt.annotate_rows(**{
            case_het_gt_ab: hl.agg.filter(case_filter,
                                          hl.cond((mt[case_het_gt_count] > 0) & hl.is_defined(mt[case_het_gt_count]),
                                                  hl.agg.count_where(passing_het_gts) / mt[case_het_gt_count],
                                                  hl.null(hl.tfloat)))})

        # Get control het GT count
        cont_filter = (mt[args.pheno_col] == False) & (hl.is_defined(mt[args.pheno_col])) &\
                      (mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) == 0) & \
                      hl.is_defined(mt.pop_outlier_sample) & hl.is_defined(mt.failing_samples_qc)
        het_gt_filters = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt[failing_dp_name]) == 0) & \
                         hl.is_defined(mt[failing_dp_name])

        mt = mt.annotate_rows(**{cont_het_gt_count: hl.agg.filter(cont_filter, hl.agg.count_where(het_gt_filters))})

        # Get control het GT fraction in allelic balance
        cont_filter = (mt[args.pheno_col] == False) & (hl.is_defined(mt[args.pheno_col])) &\
                      (mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) == 0) & \
                      hl.is_defined(mt.pop_outlier_sample) & hl.is_defined(mt.failing_samples_qc)
        passing_het_gts = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt[failing_dp_name]) == 0) & \
                          ~(mt[failing_ab_name].contains("failing_het_ab")) & hl.is_defined(mt[failing_dp_name])

        mt = mt.annotate_rows(**{
            cont_het_gt_ab: hl.agg.filter(cont_filter,
                                          hl.cond((mt[cont_het_gt_count] > 0) & hl.is_defined(mt[cont_het_gt_count]),
                                                  hl.agg.count_where(passing_het_gts) / mt[cont_het_gt_count],
                                                  hl.null(hl.tfloat)))})

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

    return mt


def find_failing_variants(mt, args, mode):
    """
    Function to find variants failing on QC measures, which can be run in 'low_pass' or 'final' mode, with varying
    filters given by args depending on the mode.

    :param mt: matrix table to filter
    :param args: arguments object with thresholds
    :return: returns matrix table with entries and variants annotated with info on whether they fail QC metrics.
    """
    ###########################################################
    # Define specific things for low pass vs final variant QC #
    ###########################################################
    if mode == 'low_pass':
        # Check that multiallelic variants have been split
        try:
            test = hl.is_defined(mt.row.was_split)
        except Exception as e:
            logging.info("Split multi-allelics before running!")
            logging.info(e)
            return

        # Define variables
        p_hwe = args.low_pass_p_hwe
        call_rate = args.low_pass_min_call_rate
        annotation_name = "variant_qc_thresholds_lowpass"
        sex_aware_call_rate = "False"
        varqc_name = 'prefilter_variant_qc'

    elif mode == "final":
        # Check that samples QC has been run
        try:
            test = hl.is_defined(mt.failing_samples_qc)
            test = hl.is_defined(mt.pop_outlier_sample)
        except Exception as e:
            logging.info('failing_samples_qc and pop_outlier_sample not defined! Run samples QC before invoking this.')
            logging.info(e)

        # Define variables
        p_hwe = args.final_p_hwe
        call_rate = args.final_min_call_rate
        annotation_name = "variant_qc_thresholds_final"
        sex_aware_call_rate = "True"
        varqc_name = 'final_no_failing_samples_varqc'

    else:
        logging.error("Error! modes allowed are low_pass and final variant QC")
        return

    failing_name = mode + "_failing_variant_qc"
    #####################################
    # Annotate QC thresholds to globals #
    #####################################
    if args.pheno_col is not None:
        hwe_tag = "controls samples only"
    else:
        hwe_tag = "all samples"

    threshold_dict = {'pval_hwe': str(p_hwe), 'hwe_excluded_in': str(hwe_tag),
                      'call_rate': str(call_rate),
                      'sex_aware_call_rate': sex_aware_call_rate, 'snp_qd': str(args.snp_qd),
                      'indel_qd': str(args.indel_qd),
                      'het_max_ref_reads_thresh': str(args.max_het_ref_reads),
                      'het_min_ref_reads_thresh': str(args.min_het_ref_reads)}

    mt = mt.annotate_globals(**{annotation_name: threshold_dict})

    ############################################################
    # Find failing genotypes and do allelic balance annotation #
    ############################################################
    # Find genotypes failing on depth, quality, and allelic balance
    mt = find_failing_genotypes_depth_quality(mt, args, mode)
    mt = find_failing_genotypes_ab(mt, args, mode)

    # Define het allelic balance value (percentage of GT calls for that variant *in balance*)
    mt = annotate_variant_het_ab(mt, args, mode)
    mt = mt.checkpoint(args.output_stem + mode + '_filtervartemp1_deleteme.mt', overwrite=True)

    ###########################################################
    # Filter matrix table to only passing genotypes + samples #
    ###########################################################
    # We have to actually filter out failing genotypes and samples because hl.variant_qc has no filter option
    # Filter out failing genotypes for calculating variant QC statistics
    logging.info("Filtering out genotypes failing on depth, quality, and allelic balance metrics for QC "
                 "metric calculation. Variant-level % het GTs in allelic balance has already been annotated,"
                 "after removing genotypes failing on depth and quality, so those variant annotations remain.")
    gtcount = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    mt_filt = sq.filter_failing(mt, args, mode, variants=False)

    filt_gt_count = mt_filt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt_filt.GT)))
    logging.info(f"Number of genotypes on which variant QC measures calculated: {filt_gt_count}")
    logging.info(f"({100 - round((gtcount - filt_gt_count) / gtcount, 2)}% of all genotypes)")

    # Run variant QC on filtered matrix table
    mt_filt = hl.variant_qc(mt_filt, name=varqc_name)

    ######################################################
    # Set HWE calculation in either all or only controls #
    ######################################################
    if args.pheno_col is not None:
        case_filters = (mt_filt[args.pheno_col] == False) & hl.is_defined(mt_filt[args.pheno_col])
        mt_filt = mt_filt.annotate_rows(
            hwe_ctrls_only=hl.agg.filter(case_filters, hl.agg.hardy_weinberg_test(mt_filt.GT)))
        logging.info('Calculated HWE with controls only (ignoring samples with missing phenotype information)')
    else:
        logging.info('Phenotype column name not provided, calculating HWE for variants using all samples.')

    mt_filt = mt_filt.checkpoint(args.output_stem + mode + '_filtervartemp3_deleteme.mt', overwrite=True)

    ###################
    # Filter variants #
    ###################
    # Initializes a row annotation as an empty array, and adds str elements if variant failing a test
    # hl.cond takes the condition, and returns the array plus a tag if the condition is satisfied,
    # or else it returns the array as it was.
    filter_dict = {}

    # Initalize row annotation with variant fail tags #
    mt_filt = mt_filt.annotate_rows(**{failing_name: hl.empty_array(hl.tstr)})
    total_variants = mt_filt.count_rows()

    ########################
    # Failing GSQR filters #
    ########################
    # Define filter condition
    logging.info("Marking variants with GQSR filters not empty")
    gqse_filter = hl.cond((hl.len(mt_filt.filters) != 0) & hl.is_defined(mt_filt.filters),
                          mt_filt[failing_name].append("failing_GQSR_filters"),
                          mt_filt[failing_name])
    # Annotate rows
    mt_filt = mt_filt.annotate_rows(**{failing_name: gqse_filter})
    # Get number of failing variants
    failing_gqsr = mt_filt.aggregate_rows(hl.agg.count_where(mt_filt[failing_name].contains("failing_GQSR_filters")))
    failing_gqsr_perc = round(failing_gqsr / total_variants * 100, 2)
    filter_dict["failing_GSQR"] = failing_gqsr
    # Report if variant annotation undefined
    gqsr_defined = mt_filt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt_filt.filters)))
    if gqsr_defined == 0:
        logging.error(f"Note! mt.filters annotation undefined for all variants! Not filtering variants on this measure.")
    else:
        gqse_miss = hl.cond(~hl.is_defined(mt_filt.filters), mt_filt[failing_name].append("missing_GQSR_filters"),
                            mt_filt[failing_name])

        mt_filt = mt_filt.annotate_rows(**{failing_name: gqse_miss})

        missing_gqsr = mt_filt.aggregate_rows(hl.agg.count_where(mt_filt[failing_name].contains("missing_GQSR_filters")))
        filter_dict['missing_GQSR'] = missing_gqsr
        missing_gqsr_perc = round(missing_gqsr / total_variants * 100, 2)

    #####################################################################
    # Find variants that are snps or indels with QD less than threshold #
    #####################################################################
    # Define filter condition
    snps = hl.is_snp(mt_filt.row.alleles[0], mt_filt.row.alleles[1]) & (mt_filt.row.info.QD < args.snp_qd) & \
           hl.is_defined(mt_filt.row.info.QD)
    indels = (~hl.is_snp(mt_filt.row.alleles[0], mt_filt.row.alleles[1])) & (mt_filt.row.info.QD < args.indel_qd) & \
             hl.is_defined(mt_filt.row.info.QD)
    qd_filter = hl.cond(snps | indels, mt_filt[failing_name].append("failing_QD"), mt_filt[failing_name])
    # Annotate rows
    mt_filt = mt_filt.annotate_rows(**{failing_name: qd_filter})
    # Get number of failing variants
    failing_qd = mt_filt.aggregate_rows(hl.agg.count_where(mt_filt[failing_name].contains("failing_QD")))
    failing_qd_perc = round(failing_qd / total_variants * 100, 2)
    filter_dict["failing_QD"] = failing_qd
    # Report if variant annotation undefined
    qd_defined = mt_filt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt_filt.row.info.QD)))
    if qd_defined < total_variants:
        qd_miss = hl.cond(~(hl.is_defined(mt_filt.row.info.QD)), mt_filt[failing_name].append("missing_QD"),
                          mt_filt[failing_name])

        mt_filt = mt_filt.annotate_rows(**{failing_name: qd_miss})

        missing_qd = mt_filt.aggregate_rows(hl.agg.count_where(mt_filt[failing_name].contains("missing_QD")))
        filter_dict['missing_QD'] = missing_qd
        missing_qd_perc = round(missing_qd / total_variants * 100, 2)

    ###################################################################
    # Find variants with >20% of het genotypes out of allelic balance #
    ###################################################################
    # Define filter condition
    logging.info("Finding variants failing on het GT allelic balance")
    failing_het_ab_name = mode + '_frac_het_gts_in_ab'
    ab_filter = hl.cond(hl.is_defined(mt_filt[failing_het_ab_name]) &
                        (mt_filt[failing_het_ab_name] < args.ab_allowed_dev_het),
                        mt_filt[failing_name].append("failing_het_ab"), mt_filt[failing_name])
    # Annotate rows
    mt_filt = mt_filt.annotate_rows(**{failing_name: ab_filter})
    # Get number of failing variants
    failing_ab = mt_filt.aggregate_rows(hl.agg.count_where(mt_filt[failing_name].contains("failing_het_ab")))
    failing_ab_perc = round(failing_ab / total_variants * 100, 2)
    filter_dict["failing_het_allelic_balance"] = failing_ab
    # Report if variant annotation undefined
    ab_defined = mt_filt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt_filt[failing_het_ab_name])))
    if ab_defined < total_variants:
        # Some variants might not have het genotypes, so we can skip filtering the variants with ab not defined.
        logging.error(f"Note! mt.{failing_het_ab_name} annotation defined for only {ab_defined} variants! "
                      f"Variants missing this annotation not filtered on this measure.")

    ####################
    # Failing HWE test #
    ####################
    # Define filter condition
    logging.info(f"Finding variants failing HWE test in {hwe_tag}")
    if args.pheno_col is not None:
        hwe_cond = (mt_filt.row.hwe_ctrls_only.p_value < p_hwe) & hl.is_defined(mt_filt.hwe_ctrls_only.p_value)
        hwe_defined = mt_filt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt_filt.hwe_ctrls_only.p_value)))
    else:
        hwe_cond = (mt_filt.row[varqc_name].p_value_hwe < p_hwe) & hl.is_defined(mt_filt[varqc_name].p_value_hwe)
        hwe_defined = mt_filt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt_filt[varqc_name].p_value_hwe)))
    hwe_filter = hl.cond(hwe_cond, mt_filt[failing_name].append("failing_hwe"), mt_filt[failing_name])
    # Annotate rows
    mt_filt = mt_filt.annotate_rows(**{failing_name: hwe_filter})
    # Get number of failing variants
    failing_hwe = mt_filt.aggregate_rows(hl.agg.count_where(mt_filt[failing_name].contains('failing_hwe')))
    failing_hwe_perc = round(failing_hwe / total_variants * 100, 2)
    filter_dict["failing_hwe"] = failing_hwe
    # Report if variant annotation undefined
    if hwe_defined < total_variants:
        # We should always expect defined values for this.
        logging.error(f"Note! HWE annotation defined for only {hwe_defined} variants! "
                      f"Something is wrong, check this.")

    ##############################################
    # Find variants not passing call rate filter #
    ##############################################
    # Define filter condition
    if mode == "final":
        logging.info("Finding variants not passing call rate filter, using sex-aware variant call rate ")
        call_rate_cond = (mt_filt.row.sexaware_call_rate < args.final_min_call_rate) & \
                          hl.is_defined(mt_filt.sexaware_call_rate)
        cr_defined = mt_filt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt_filt.sexaware_call_rate)))
        low_passing_vars = mt_filt.aggregate_rows(hl.agg.count_where(hl.len(mt_filt.low_pass_failing_variant_qc) == 0))
        if cr_defined < low_passing_vars:
            logging.error(f"Note! sex aware call rate annotation defined for only {cr_defined} variants, where "
                          f"there are {low_passing_vars} variants passing low pass variant QC.Something is wrong!")

    else:  # low-pass
        logging.info("Finding variants not passing call rate filter, NOT using sex-aware variant call rate ")
        call_rate_cond = (mt_filt[varqc_name].call_rate < args.low_pass_min_call_rate) & \
                          hl.is_defined(mt_filt[varqc_name].call_rate)
        cr_defined = mt_filt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt_filt[varqc_name].call_rate)))
        if cr_defined < total_variants:
            # We should always expect defined values for this.
            logging.error(f"Note! Call rate annotation defined for only {cr_defined} variants! "
                          f"Something is wrong, check this.")

    call_rate_filter = hl.cond(call_rate_cond, mt_filt[failing_name].append("failing_call_rate"), mt_filt[failing_name])

    # Annotate rows
    mt_filt = mt_filt.annotate_rows(**{failing_name: call_rate_filter})
    # Get number of failing variants
    failing_call_rate = mt_filt.aggregate_rows(hl.agg.count_where(mt_filt[failing_name].contains("failing_call_rate")))
    failing_cr_perc = round(failing_call_rate / total_variants * 100, 2)
    filter_dict["failing_call_rate"] = failing_call_rate

    ######################################################
    # Annotate failing counts to globals, report to logs #
    ######################################################
    # Report exclusion numbers to logs
    logging.info(f"\nVariants with length mt.row.filters != 0: {failing_gqsr}, {failing_gqsr_perc}%"
                 f"\nsnp QD < {args.snp_qd} or indel QD < {args.indel_qd}: {failing_qd}, {failing_qd_perc}%"
                 f"\n>{args.ab_allowed_dev_het*100}% of het genotypes out of allelic balance: "
                 f"{failing_ab}, {failing_ab_perc}%"
                 f"\np value HWE < {p_hwe} in {hwe_tag}: {failing_hwe}, {failing_hwe_perc}%"
                 f"\ncall rate < {getattr(args, mode + '_min_call_rate')}: {failing_call_rate}, {failing_cr_perc}%")

    passing = mt_filt.aggregate_rows(hl.agg.count_where(hl.len(mt_filt[failing_name]) == 0))

    logging.info(f"Variants passing QC: {passing}, {round(passing / total_variants * 100, 2)}%")
    # Add annotations back to main matrix table (not filtered), return that
    mt = mt.annotate_globals(**{mode + "_filtered_counts": filter_dict})
    mt = mt.annotate_rows(**{varqc_name: mt_filt.index_rows(mt.row_key)[varqc_name]})
    mt = mt.annotate_rows(**{failing_name: mt_filt.index_rows(mt.row_key)[failing_name]})

    if args.pheno_col is not None:
        mt = mt.annotate_rows(hwe_ctrls_only=mt_filt.index_rows(mt.row_key).hwe_ctrls_only)

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

    ##############################################################
    # Annotate variants failing on case-specific allelic balance #
    ##############################################################
    case_cond = hl.cond(hl.is_defined(mt.final_case_frac_het_gts_in_ab) &
                        (mt.final_case_frac_het_gts_in_ab < args.ab_allowed_dev_het),
                        mt.failing_pheno_varqc.append("failing_case_het_ab"),
                        mt.failing_pheno_varqc)
    mt = mt.annotate_rows(failing_pheno_varqc=case_cond)
    failing_ab_case = mt.aggregate_rows(hl.agg.count_where(mt.failing_pheno_varqc.contains("failing_case_het_ab")))

    control_cond = hl.cond(hl.is_defined(mt.final_control_frac_het_gts_in_ab) &
                           (mt.final_control_frac_het_gts_in_ab < args.ab_allowed_dev_het),
                           mt.failing_pheno_varqc.append("failing_control_het_ab"),
                           mt.failing_pheno_varqc)
    mt = mt.annotate_rows(failing_pheno_varqc=control_cond)
    failing_ab_cont = mt.aggregate_rows(hl.agg.count_where(mt.failing_pheno_varqc.contains("failing_control_het_ab")))

    ################################################################
    # Find variants failing on case-specific call rate (sex aware) #
    ################################################################
    case_cr_cond = hl.cond(hl.is_defined(mt.sexaware_case_call_rate) &
                           (mt.sexaware_case_call_rate < args.pheno_call_rate),
                           mt.failing_pheno_varqc.append("failing_case_call_rate"),
                           mt.failing_pheno_varqc)
    mt = mt.annotate_rows(failing_pheno_varqc=case_cr_cond)
    failing_cr_case = mt.aggregate_rows(hl.agg.count_where(mt.failing_pheno_varqc.contains("failing_case_call_rate")))

    cont_cr_cond = hl.cond(hl.is_defined(mt.sexaware_cont_call_rate) &
                           (mt.sexaware_cont_call_rate < args.pheno_call_rate),
                           mt.failing_pheno_varqc.append("failing_cont_call_rate"),
                           mt.failing_pheno_varqc)
    mt = mt.annotate_rows(failing_pheno_varqc=cont_cr_cond)
    failing_cr_cont = mt.aggregate_rows(hl.agg.count_where(mt.failing_pheno_varqc.contains("failing_cont_call_rate")))

    mt = mt.annotate_globals(case_control_callrate_threshold=args.pheno_call_rate)

    logging.info(f"Number of variants failing on case ab: {failing_ab_case}")
    logging.info(f"Number of variants failing on control ab: {failing_ab_cont}")
    logging.info(f"Number of variants failing on case call rate: {failing_cr_case}")
    logging.info(f"Number of variants failing on control call rate: {failing_cr_cont}")

    failing_any = mt.aggregate_rows(hl.agg.count_where(hl.len(mt.failing_pheno_varqc) != 0))
    logging.info(f"Number of variants failing on any phenotype-specific measure: {failing_any}")

    return mt
