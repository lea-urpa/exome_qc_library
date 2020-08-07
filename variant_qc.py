"""
Functions for variant QC for exome sequencing data with Hail.

Author: Lea Urpa, August 2020
"""
import logging
import hail as hl


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
    mt = mt.annotate_entries(**{failing_name: hl.empty_array(hl.tstr)})

    # Filter entries for depth
    depth_cond = hl.cond(mt.DP < args.min_dp, mt.failing_depth_quality.append("failing_DP"), mt.failing_depth_quality)
    mt = mt.annotate_entries(**{failing_name: depth_cond})
    failing_depth = mt.aggregate_entries(hl.agg.count_where(mt.failing_depth_quality.contains("failing_DP")))
    failing_depth_perc = round(failing_depth / gt_total*100, 2)

    # Filter GQ or PL for different genotypes
    gq_cond = hl.cond((mt.GT.is_het() | mt.GT.is_hom_var()) & (mt.PL[0] < args.min_gq),
                      mt.failing_depth_quality.append("failing_PL"), mt.failing_depth_quality)
    mt = mt.annotate_entries(**{failing_name: gq_cond})

    pl_cond = hl.cond(mt.GT.is_hom_ref() & (mt.GQ < args.min_gq), mt.failing_depth_quality.append("failing_GQ"),
                      mt.failing_depth_quality)
    mt = mt.annotate_entries(**{failing_name: pl_cond})

    failing_gq = mt.aggregate_entries(hl.agg.count_where(mt.failing_depth_quality.contains("failing_PL") |
                                                         mt.failing_depth_quality.contains("failing_GQ")))
    failing_gq_perc = round(failing_gq / gt_total * 100, 2)

    logging.info(f"Number of GTs with depth < {args.min_dp}: {failing_depth} ({failing_depth_perc}%)")
    logging.info(f"Number of GTs with GQ/PL < {args.min_gq}: {failing_gq} ({failing_gq_perc}%)")

    passing_gts = mt.aggregate_entries(hl.agg.count_where(hl.len(mt.failing_depth_quality) == 0))
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

    mt = mt.annotate_globals(genotype_qc_thresholds_ab={'max_het_ref_reads': args.max_het_ref_reads,
                                                        'min_het_ref_reads': args.min_het_ref_reads,
                                                        'min_hom_ref_ref_reads': args.min_hom_ref_ref_reads,
                                                        'max_hom_alt_ref_reads': args.max_hom_alt_ref_reads})

    ##################################################################
    # Check that depth and gq filters have been applied, then filter #
    ##################################################################
    try:
        test = hl.is_defined(mt.failing_depth_quality)
    except Exception as e:
        logging.error("Error! Run find_failing_genotypes_depth_quality() on this matrix table before running "
                      "find_failing_genotypes_ab().")
        logging.error(e)
        exit()

    ################################
    # Get starting genotypes count #
    ################################
    mt = mt.annotate_entries(failing_ab=hl.empty_array(hl.tstr))
    gthet = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_het()))
    gthomvar = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_hom_var()))
    gthomref = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_hom_ref()))

    ###########################################
    # Filter het genotypes on allelic balance #
    ###########################################
    het_ab_cond = ((mt.AD[0] / hl.sum(mt.AD)) < args.min_het_ref_reads) | \
                  ((mt.AD[0] / hl.sum(mt.AD)) > args.max_het_ref_reads)
    ab_cond_het = hl.cond(mt.GT.is_het() & (hl.len(mt.failing_depth_quality) == 0) & het_ab_cond,
                          mt.failing_ab.append("failing_het_ab"), mt.failing_ab)

    mt = mt.annotate_entries(failing_ab=ab_cond_het)

    hets_failing_ab = mt.aggregate_entries(hl.agg.count_where(mt.failing_ab.contains("failing_het_ab")))
    het_failing_ab_perc = round(hets_failing_ab / gthet * 100, 2)

    ###############################################
    # Filter hom ref genotypes on allelic balance #
    ###############################################
    hom_ab_cond = (mt.AD[0] / hl.sum(mt.AD)) < args.min_hom_ref_ref_reads
    ab_cond_homref = hl.cond(mt.GT.is_hom_ref() & hom_ab_cond, mt.failing_ab.append("failing_homref_ab"), mt.failing_ab)

    mt = mt.annotate_entries(failing_ab=ab_cond_homref)

    homref_failing_ab = mt.aggregate_entries(hl.agg.count_where(mt.failing_ab.contains("failing_homref_ab")))
    homref_failing_ab_perc = round(homref_failing_ab / gthomref * 100, 2)

    ###############################################
    # Filter hom var genotypes on allelic balance #
    ###############################################
    homalt_ab_cond = (mt.AD[0] / hl.sum(mt.AD)) > args.max_hom_alt_ref_reads
    ab_cond_homalt = hl.cond(mt.GT.is_hom_var() & homalt_ab_cond, mt.failing_ab.append("failing_homalt_ab"),
                             mt.failing_ab)

    mt = mt.annotate_entries(failing_ab=ab_cond_homalt)

    homalt_failing_ab = mt.aggregate_entries(hl.agg.count_where(mt.failing_ab.contains("failing_homalt_ab")))
    homalt_failing_ab_perc = round(homalt_failing_ab / gthomvar * 100, 2)

    ###################################################
    # Report number and percent of genotypes excluded #
    ###################################################
    logging.info("Genotypes failing ab, after excluding low DP and low GQ genotypes:")
    logging.info(f"Number of het GTs failing ab: {hets_failing_ab}({het_failing_ab_perc}%)")
    logging.info(f"Number of hom ref GTs excluded: {homref_failing_ab} ({homref_failing_ab_perc}%)")
    logging.info(f"Number of hom alt GTs excluded: {homalt_failing_ab} ({homalt_failing_ab_perc})")

    mt = mt.annotate_globals(genotype_qc_exclusions=
                             {'het_excluded_ct_percent': [hets_failing_ab, het_failing_ab_perc],
                              'hom_ref_excluded_ct_percent': [homref_failing_ab, homref_failing_ab_perc],
                              'hom_var_excluded_ct_percent': [homalt_failing_ab, homalt_failing_ab_perc]})

    return mt


def annotate_variant_het_ab(mt, args, prefix=""):
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
                 "population outliers and samples failing samples QC (if run already)")
    het_gt_cnt = prefix + 'het_gt_count'
    het_ab = prefix + 'frac_het_gts_in_ab'

    # Check that genotype filter annotations have been added
    try:
        hl.is_defined(mt.failing_ab)
        hl.is_defined(mt.failing_depth_quality)
    except Exception as e:
        logging.error("Error! find_failing_genotypes_ab() and find_failing_genotypes_depth_quality() should be run "
                      "before annotating variant's % of hets in allelic balance, to filter out bad genotypes.")
        logging.error(e)
        exit()

    # Check if samples QC has been run yet, if not add empty annotations
    try:
        hl.is_defined(mt.pop_outlier_sample)
    except Exception as e:
        logging.info("Detected samples QC not run, adding False pop_outlier_sample annotation to columns.")
        logging.info(e)
        mt = mt.annotate_cols(pop_outlier_sample=False)

    try:
        hl.is_defined(mt.failing_samples_qc)
    except Exception as e:
        logging.info("Detected samples QC not run, adding empty array failing_samples_qc annotation to columns.")
        logging.info(e)
        mt = mt.annotate_cols(failing_samples_qc=hl.empty_array(hl.tstr))

    ################################################################
    # Annotate rows with het GT count, then het GT allelic balance #
    ################################################################
    # Get count of het genotypes per variant
    case_filter = (mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) == 0)
    het_gt_filters = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt.failing_depth_quality) == 0)

    mt = mt.annotate_rows(**{het_gt_cnt: hl.agg.filter(case_filter, hl.float(hl.agg.count_where(het_gt_filters)))})

    # Get percent of het gts in AB if het gt count != 0
    case_filter = (mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) == 0)
    passing_het_gts = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt.failing_depth_quality) == 0) & \
                      ~(mt.failing_ab.contains("failing_het_ab"))

    # Counting het GTs skips variants where there are no het GTs (which leads to div by 0 error)
    mt = mt.annotate_rows(**{het_ab:
             hl.agg.filter(case_filter,
                           hl.cond(mt[het_gt_cnt] > 0,
                                   hl.float(hl.agg.count_where(passing_het_gts)) / mt[het_gt_cnt],
                                   hl.null(hl.tfloat)))})

    ################################################################################
    # Annotate het GT and het GT allelic balance for cases and controls separately #
    ################################################################################
    if args.pheno_col is not None:
        case_het_gt_count = prefix + 'case_het_gt_count'
        case_het_gt_ab = prefix + 'case_frac_het_gts_in_ab'
        cont_het_gt_count = prefix + 'control_het_gt_count'
        cont_het_gt_ab = prefix + 'control_frac_het_gts_in_ab'

        # Count case het GTs
        case_filter = (mt[args.pheno_col] == True) & (mt.pop_outlier_sample == False) & \
                      (hl.len(mt.failing_samples_qc) == 0)
        het_gt_filters = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt.failing_depth_quality) == 0)
        mt = mt.annotate_rows(**{case_het_gt_count: hl.agg.filter(case_filter,
                                                                  hl.agg.count_where(het_gt_filters))})

        # Get case het GT fraction in allelic balance
        passing_het_gts = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt.failing_depth_quality) == 0) & \
                          ~(mt.failing_ab.contains("failing_het_ab"))

        case_filter = (mt[args.pheno_col] == True) & (mt.pop_outlier_sample == False) & \
                      (hl.len(mt.failing_samples_qc > 0))

        mt = mt.annotate_rows(**{
            case_het_gt_ab: hl.agg.filter(case_filter,
                                          hl.cond(mt[case_het_gt_count] > 0,
                                                  hl.float(hl.agg.count_where(passing_het_gts)) / mt[case_het_gt_count],
                                                  hl.null(hl.tfloat)))})

        # Get control het GT count
        cont_filter = (mt[args.pheno_col] == False) & (mt.pop_outlier_sample == False) & \
                      (hl.len(mt.failing_samples_qc) == 0)
        het_gt_filters = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt.failing_depth_quality) == 0)

        mt = mt.annotate_rows(**{cont_het_gt_count: hl.agg.filter(cont_filter,
                                                                  hl.float(hl.agg.count_where(het_gt_filters)))})

        # Get control het GT fraction in allelic balance
        cont_filter = (mt[args.pheno_col] == False) & (mt.pop_outlier_sample == False) & \
                      (hl.len(mt.failing_samples_qc) == 0)
        passing_het_gts = mt.GT.is_het() & hl.is_defined(mt.GT) & (hl.len(mt.failing_depth_quality) == 0) & \
                          ~(mt.failing_ab.contains("failing_het_ab"))

        mt = mt.annotate_rows(**{
            cont_het_gt_ab: hl.agg.filter(cont_filter,
                                          hl.cond(mt[cont_het_gt_count] > 0,
                                                  hl.float(hl.agg.count_where(passing_het_gts)) / mt[cont_het_gt_count],
                                                  hl.null(hl.tbool)))})

    return mt


def find_failing_variants(mt, args, mode):
    """
    Function to find variants failing on QC measures, which can be run in 'low_pass' or 'final' mode, with varying
    filters given by args depending on the mode.

    :param mt: matrix table to filter
    :param args: arguments object with thresholds
    :return: returns matrix table with variants filtered out.
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
        failing_name = args.lowpass_fail_name

    elif mode == "final":
        # Check that samples QC has been run
        try:
            test = hl.is_defined(mt.failing_qc_samples)
            test = hl.is_defined(mt.pop_outlier_samples)
        except Exception as e:
            logging.info('failing_qc_samples and pop_outlier_samples not defined! Run samples QC before invoking this.')
            logging.info(e)

        # Define variables
        p_hwe = args.final_p_hwe
        call_rate = args.final_min_call_rate
        annotation_name = "variant_qc_thresholds_final"
        sex_aware_call_rate = "True"
        varqc_name = 'final_no_failing_samples_varqc'
        failing_name = args.final_fail_name

    else:
        logging.error("Error! modes allowed are low_pass and final variant QC")
        return

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
    mt = annotate_variant_het_ab(mt, args, prefix=mode)
    mt = mt.checkpoint(args.output_stem + mode + '_filtervartemp1_deleteme.mt', overwrite=True)

    ###########################################################
    # Filter matrix table to only passing genotypes + samples #
    ###########################################################
    # We have to actually filter out failing genotypes and samples because hl.variant_qc has no filter option
    # Filter out failing genotypes for calculating variant QC statistics
    logging.info("Filtering out genotypes failing on depth, quality, and allelic balance metrics for initial QC "
                 "metric calculation.")
    gtcount = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    mt_filt = mt.filter_entries(hl.len(mt.failing_depth_quality) == 0, keep=True)
    mt_filt = mt_filt.filter_entries(hl.len(mt.failing_ab) == 0, keep=True)

    mt_filt = mt_filt.checkpoint(args.output_stem + mode + "filter1_deleteme.mt")
    filt_gt_count = mt_filt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    logging.info(f"Number of genotypes on which initial variant QC measures calculated: {filt_gt_count}")
    logging.info(f"({100 - round(gtcount - filt_gt_count, 2)}% of all genotypes)")

    # Filter out samples failing QC for calculating variant QC statistics (empty if lowpass)
    logging.info("Filtering out samples failing samples QC. Annotations are empty if samples QC not run yet.")
    samples_count = mt.count_cols()
    mt_filt = mt_filt.filter_cols((mt.pop_outlier_sample == False) & (hl.len(mt.failing_samples_qc) == 0))
    mt_filt = mt_filt.checkpoint(args.output_stem + mode + "filter2_deleteme.mt")
    filt_s_count = mt_filt.count_cols()
    logging.info(f"Initial samples count {samples_count}. Sample count after filtering: {filt_s_count}")

    # Run variant QC on filtered matrix table
    mt_filt = hl.variant_qc(mt_filt, name=varqc_name)

    if mode == 'final':
        gt_filters = hl.is_defined(mt_filt.GT) & mt_filt.GT.is_het()
        mt_filt = mt_filt.annotate_rows(final_het_callrate=hl.agg.count_where(gt_filters) / mt_filt.final_het_gt_count)

        # Sanity check: is the het callrate ever higher than het allelic balance?
        logging.info('Any variants where het callrate is higher than het allelic balance? Should be zero.')
        logging.info(mt_filt.aggregate_rows(hl.agg.count_where(
            mt_filt.final_het_callrate > mt_filt.final_frac_het_gts_in_ab)))

    ######################################################
    # Set HWE calculation in either all or only controls #
    ######################################################
    #TODO continue adapting the code below to mt_filt instead of mt, collecting the annotations that need to be added
    # back at the end. Also, how the case filters in the rest of the code deal with missing pheno col data.

    if args.pheno_col is not None:
        case_filters = (mt_filt[args.pheno_col] == False) | hl.is_missing(mt_filt.col[args.pheno_col])
        mt_filt = mt_filt.annotate_rows(
            hwe_ctrls_only=hl.agg.filter(case_filters, hl.agg.hardy_weinberg_test(mt_filt.GT)))
        logging.info('Annotated with hwe_ctrls_only successfully')
    else:
        logging.info('Phenotype column name not provided, calculating HWE for variants using all samples.')

    mt = mt.checkpoint(args.output_stem + '_filtervartemp3_deleteme.mt', overwrite=True)

    ###################
    # Filter variants #
    ###################
    # Initializes a row annotation as an empty array, and adds str elements if variant failing a test
    # hl.cond takes the condition, and returns the array plus a tag if the condition is satisfied,
    # or else it returns the array as it was.
    filter_dict = {}

    # Initalize row annotation with variant fail tags #
    mt = mt.annotate_rows(**{failing_name: hl.empty_array(hl.tstr)})
    total_variants = mt.count_rows()

    # Failing GSQR filters #
    logging.info("Marking variants with GQSR filters not empty")
    gqse_filter = hl.cond(hl.len(mt.filters != 0), mt[failing_name].append("failing_GQSR_filters"), mt[failing_name])
    mt = mt.annotate_rows(**{failing_name: gqse_filter})
    failing_gqsr = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains["failing_GQSR_filters"]))
    failing_gqsr_perc = round(failing_gqsr / total_variants * 100, 2)
    filter_dict["failing_GSQR"] = failing_gqsr

    # Find variants that are snps or indels with QD less than threshold #
    snps = hl.is_snp(mt.row.alleles[0], mt.row.alleles[1]) & (mt.row.info.QD < args.snp_qd)
    indels = (~hl.is_snp(mt.row.alleles[0], mt.row.alleles[1])) & (mt.row.info.QD < args.indel_qd)
    qd_filter = hl.cond(snps | indels, mt[failing_name].append("failing_QD"), mt[failing_name])
    mt = mt.annotate_rows(**{failing_name: qd_filter})
    failing_qd = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains("failing_QD")))
    failing_qd_perc = round(failing_qd / total_variants * 100, 2)
    filter_dict["failing_QD"] = failing_qd

    # Find variants with >20% of het genotypes out of allelic balance #
    logging.info("Finding variants failing on het GT allelic balance")
    ab_filter = hl.cond(hl.is_defined(mt.final_het_ab_all) & (mt.final_het_ab_all < args.final_ab_allowed_dev_het),
                        mt[failing_name].append("failing_het_ab"), mt[failing_name])
    mt = mt.annotate_rows(**{failing_name: ab_filter})
    failing_ab = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains("failing_het_ab")))
    failing_ab_perc = round(failing_ab / total_variants * 100, 2)
    filter_dict["failing_het_allelic_balance"] = failing_ab

    # Failing HWE test #
    logging.info(f"Finding variants failing HWE test in {hwe_tag}")
    if args.pheno_col is not None:
        hwe_cond = mt.row.hwe_ctrls_only.p_value < p_hwe
    else:
        hwe_cond = mt.row[varqc_name].p_value_hwe < p_hwe  # var QC name given at top of function

    hwe_filter = hl.cond(hwe_cond, mt[failing_name].append("failing_hwe"), mt[failing_name])
    mt = mt.annotate_rows(**{failing_name: hwe_filter})
    failing_hwe = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains['failing_hwe']))
    failing_hwe_perc = round(failing_hwe / total_variants * 100, 2)
    filter_dict["failing_hwe"] = failing_hwe

    # Find variants not passing call rate filter #
    if mode == "final":
        logging.info("Finding variants not passing call rate filter, using sex-aware call rate ")
        call_rate_cond = mt.row.sexaware_call_rate < args.min_call_rate

    else:  # low-pass
        logging.info("Finding variants not passing call rate filter, NOT using sex-aware call rate ")
        call_rate_cond = mt[varqc_name].call_rate < args.min_call_rate

    call_rate_filter = hl.cond(call_rate_cond, mt[failing_name].append("failing_call_rate"), mt[failing_name])
    mt = mt.annotate_rows(**{failing_name: call_rate_filter})
    failing_call_rate = mt.aggregate_rows(hl.agg.count_where(mt[failing_name].contains("failing_call_rate")))
    failing_cr_perc = round(failing_call_rate / total_variants * 100, 2)
    filter_dict["failing_call_rate"] = failing_call_rate

    ######################################################
    # Annotate failing counts to globals, report to logs #
    ######################################################
    # Report exclusion numbers to logs
    logging.info(f"\nVariants with length mt.row.filters != 0: {failing_gqsr}, {failing_gqsr_perc}%"
                 f"\nsnp QD < {args.snp_qd} or indel QD < {args.indel_qd}: {failing_qd}, {failing_qd_perc}%"
                 f"\n>{args.final_ab_allowed_dev_het*100}% of het genotypes out of allelic balance: "
                 f"{failing_ab}, {failing_ab_perc}%"
                 f"\np value HWE < {p_hwe} in {hwe_tag}: {failing_hwe}, {failing_hwe_perc}%"
                 f"\ncall rate < {args.min_call_rate}: {failing_call_rate}, {failing_cr_perc}%")

    failing_any = mt.aggregate_rows(hl.agg.count_where(hl.len(mt[failing_name]) != 0))

    logging.info(f"Variants failing QC: {failing_any}, {round(failing_any / total_variants * 100, 2)}%")

    # Add annotations back to main matrix table (not filtered), return that
    mt = mt.annotate_rows(**{varqc_name: mt_filt.index_rows(mt.row_key)[varqc_name]})

    if mode == 'final':
        mt = mt.annotate_rows(final_het_callrate=mt_filt.index_rows(mt.row_key).final_het_callrate)




    if mode == "low_pass":
        mt = mt.annotate_globals(variant_qc_fails_lowpass=filter_dict)
        return mt

    else:  # final
        mt_orig = mt_orig.annotate_globals(variant_qc_fails_final=filter_dict)

        # Annotate original matrix table with variant failing col annotation, and allelic balance annotations
        # calculated without population outliers and samples failing QC

        mt_orig = mt_orig.annotate_rows(**{failing_name: mt.index_rows(mt_orig.row_key)[failing_name]})

        mt_orig = mt_orig.annotate_rows(**{varqc_name: mt.index_rows(mt_orig.row_key)[varqc_name]})

        mt_orig = mt_orig.annotate_rows(final_het_gt_count=mt.index_rows(mt_orig.row_key).final_het_gt_count,
                                        final_het_ab_all=mt.index_rows(mt_orig.row_key).final_het_ab_all,
                                        final_het_gt_count_cases=mt.index_rows(
                                            mt_orig.row_key).final_het_gt_count_cases,
                                        final_het_ab_cases=mt.index_rows(mt_orig.row_key).final_het_ab_cases,
                                        final_het_gt_count_cont=mt.index_rows(
                                            mt_orig.row_key).final_het_gt_count_cont,
                                        final_het_ab_cont=mt.index_rows(mt_orig.row_key).final_het_ab_cont,
                                        final_het_callrate=mt.index_rows(mt_orig.row_key).final_het_callrate)

        if args.pheno_col is not None:
            mt_orig = mt_orig.annotate_rows(hwe_ctrls_only=mt.index_rows(mt_orig.row_key).hwe_ctrls_only)

        return mt_orig


