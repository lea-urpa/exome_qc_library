"""
Functions for variant QC for exome sequencing data with Hail.

Author: Lea Urpa, August 2020
"""
import logging
import hail as hl


def filter_genotypes_depth_quality(mt, args):
    """
    Filter genotypes for initial low-pass variant QC in qc pipeline. Filters out low DP and poor GQ genotypes, and
    unfilters the entries (sets the filtered entries to missing) afterward.
    :param mt: matrix table to filter
    :return: returns filtered matrix table
    """
    # Log thresholds + report
    logging.info(f"Filtering out genotypes on min_dp <= {args.min_dp}, or min GQ <= {args.min_gq}")
    mt = mt.annotate_globals(genotype_qc_thresh_dp_gq={'min_dp': args.min_dp, 'min_gq': args.min_gq})

    # Get starting genotypes count
    gtstart = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))

    # Filter entries for depth
    mt = mt.filter_entries((mt.DP <= args.min_dp), keep=False)
    gtcount1 = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    depth_excluded = gtstart - gtcount1

    # Filter GQ or PL for different genotypes
    mt = mt.filter_entries((mt.GT.is_het() | mt.GT.is_hom_var()) & (mt.PL[0] < args.min_gq), keep=False)
    mt = mt.filter_entries(mt.GT.is_hom_ref() & (mt.GQ < args.min_gq), keep=False)
    gtcount2 = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    quality_excluded = gtcount1 - gtcount2

    depth_percent = round(depth_excluded/gtstart*100, 2)
    quality_percent = round(quality_excluded/gtstart*100, 2)

    logging.info(f"Number of GTs excluded for depth: {depth_excluded} ({depth_percent}%)")
    logging.info(f"Number of GTs excluded for quality: {quality_excluded} ({quality_percent}%)")

    mt = mt.annotate_globals(genotype_qc_lowpass_excl={'depth_excluded': depth_excluded,
                                                       'quality_excluded': quality_excluded})

    logging.info('Unfiltering entries (setting filtered entries to missing) after genotype QC.')
    mt = mt.unfilter_entries()

    return mt


def filter_genotypes_ab(mt, args):
    """
    Filters genotypes for final variant QC. Filters out genotypes with poor allelic balance, defined as the percentage
    of ref reads that should correspond to that genotype. For hets, percent of ref reads should be between 20-80%,
    for hom ref they should be at least 90%, and for hom alt they should be maximum 10%.
    :param mt: matrix table to annotated
    :param args: args to get thresholds from
    :return: 
    """
    ###############################
    # Log and annotate thresholds #
    ###############################
    logging.info(f"\nFiltering out genotypes on het ref reads > {args.max_het_ref_reads}"
                 f"\nor het ref reads < {args.min_het_ref_reads}"
                 f"\nor hom_ref ref reads < {args.min_hom_ref_ref_reads}"
                 f"\nor hom_alt ref reads > {args.max_hom_alt_ref_reads}")

    mt = mt.annotate_globals(genotype_qc_thresholds_ab={'max_het_ref_reads': args.max_het_ref_reads,
                                                        'min_het_ref_reads': args.min_het_ref_reads,
                                                        'min_hom_ref_ref_reads': args.min_hom_ref_ref_reads,
                                                        'max_hom_alt_ref_reads': args.max_hom_alt_ref_reads})
    ################################
    # Get starting genotypes count #
    ################################
    gtstart = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    logging.info('Starting count for number of genotypes: ' + str(gtstart))

    gthet = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_het()))
    gthomvar = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_hom_var()))
    gthomref = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_hom_ref()))

    ###########################################
    # Filter het genotypes on allelic balance #
    ###########################################
    mt = mt.filter_entries(mt.GT.is_het() & (((mt.AD[0] / hl.sum(mt.AD)) < args.min_het_ref_reads) |
                                             ((mt.AD[0] / hl.sum(mt.AD)) > args.max_het_ref_reads)), keep=False)
    gtcount1 = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    het_excluded = gtstart - gtcount1

    ###############################################
    # Filter hom ref genotypes on allelic balance #
    ###############################################
    mt = mt.filter_entries(mt.GT.is_hom_ref() & ((mt.AD[0] / hl.sum(mt.AD)) < args.min_hom_ref_ref_reads), keep=False)
    gtcount2 = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    hom_ref_excluded = gtcount1 - gtcount2

    ###############################################
    # Filter hom var genotypes on allelic balance #
    ###############################################
    mt = mt.filter_entries(mt.GT.is_hom_var() & ((mt.AD[0] / hl.sum(mt.AD)) > args.max_hom_alt_ref_reads), keep=False)
    gtcount3 = mt.aggregate_entries(hl.agg.count_where(hl.is_defined(mt.GT)))
    hom_var_excluded = gtcount2 - gtcount3

    ###################################################
    # Report number and percent of genotypes excluded #
    ###################################################
    het_percent = round(het_excluded/gthet*100,3)
    hom_ref_percent = round(hom_ref_excluded/gthomref*100,3)
    hom_var_percent = round(hom_var_excluded/gthomvar*100,3)
    logging.info(f"Number of het GTs excluded: {het_excluded}({het_percent}%)")
    logging.info(f"Number of hom ref GTs excluded: {hom_ref_excluded} ({hom_ref_percent}%)")
    logging.info(f"Number of hom alt GTs excluded: {hom_var_excluded} ({hom_var_percent})")

    mt = mt.annotate_globals(genotype_qc_exclusions={'het_excluded_ct_percent': [het_excluded, het_percent],
                                                     'hom_ref_excluded_ct_percent': [hom_ref_excluded, hom_ref_percent],
                                                     'hom_var_excluded_ct_percent': [hom_var_excluded, hom_var_percent]})

    logging.info('Unfiltering entries (setting filtered entries to missing) after genotype QC.')
    mt = mt.unfilter_entries()

    return mt


def annotate_ab(mt, args, prefix="", filter_bad_samples=False):
    """
     Annotate percentage of het genotypes per variant that are in  allelic balance. That means that the genotype has ref
    reads within the 20-80% range of all reads.

    :param mt: matrix table to annotate
    :param prefix: prefix to add to variant annotation
    :param filter_bad_samples: If running in final allelic balance calculation, remove 'bad' samples not passing samples
     QC?
    :return: returns annotated matrix table
    """
    #########################################################
    # Filter out population outliers and samples failing QC #
    #########################################################
    if filter_bad_samples:
        mt_orig = mt
        mt = mt.filter_cols((mt.pop_outlier_samples == True) | (mt.failing_qc_samples > 0), keep=False)

    ################################################################
    # Annotate rows with het GT count, then het GT allelic balance #
    ################################################################
    het_gt_cnt = prefix + 'het_gt_count'
    het_ab = prefix + 'het_ab_all'

    mt = mt.annotate_rows(**{het_gt_cnt: hl.float(hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT)))})

    mt = mt.annotate_rows(**{het_ab: hl.cond(mt[het_gt_cnt] > 0,  # skips variants where there are no het GTs > error
                          hl.float(hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) &
                                   ((mt.AD[0] / hl.sum(mt.AD) > args.min_het_ref_reads) &
                                    (mt.AD[0] / hl.sum(mt.AD) < args.max_het_ref_reads)))) /
                          mt[het_gt_cnt],
                          hl.null(hl.tfloat))})

    ################################################################################
    # Annotate het GT and het GT allelic balance for cases and controls separately #
    ################################################################################
    if args.pheno_col is not None:
        case_het_gt_count = prefix + 'case_het_gt_count'
        case_het_gt_ab = prefix + 'case_het_ab'
        cont_het_gt_count = prefix + 'control_het_gt_count'
        cont_het_gt_ab = prefix + 'control_het_ab'

        mt = mt.annotate_rows(**{case_het_gt_count: hl.agg.filter(mt[args.pheno_col] == True,
                                     hl.float(hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT))))})

        mt = mt.annotate_rows(**{case_het_gt_ab:
                              hl.agg.filter(mt[args.pheno_col] == True,
                                            hl.cond(mt[case_het_gt_count] > 0,
                                    hl.float(hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) &
                                                                ((mt.AD[0] / hl.sum(mt.AD) > args.min_het_ref_reads) &
                                                                 (mt.AD[0] / hl.sum(mt.AD) < args.max_het_ref_reads)))) /
                                    mt[case_het_gt_count],
                                    hl.null(hl.tfloat)))})

        mt = mt.annotate_rows(**{cont_het_gt_count: hl.agg.filter(mt[args.pheno_col] == False,
                                    hl.float(hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT))))})

        mt = mt.annotate_rows(**{cont_het_gt_ab:
                              hl.agg.filter(mt[args.pheno_col] == False,
                                            hl.cond(mt[cont_het_gt_count] > 0,
                                   hl.float(hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) &
                                                               ((mt.AD[0] / hl.sum(mt.AD) > args.min_het_ref_reads) &
                                                                (mt.AD[0] / hl.sum(mt.AD) < args.max_het_ref_reads)))) /
                                   mt[cont_het_gt_count],
                                   hl.null(hl.tbool)))})

    ###################################################################################
    # Annotate original mt with row annotations and return, or return original matrix #
    ###################################################################################
    if filter_bad_samples:
        mt_orig = mt_orig.annotate_rows(final_het_gt_count=mt.index_rows(mt_orig.row_key).final_het_gt_count,
                                        final_het_ab_all=mt.index_rows(mt_orig.row_key).final_het_ab_all,
                                        final_het_gt_count_cases=mt.index_rows(mt_orig.row_key).final_het_gt_count_cases,
                                        final_het_ab_cases=mt.index_rows(mt_orig.row_key).final_het_ab_cases,
                                        final_het_gt_count_cont=mt.index_rows(mt_orig.row_key).final_het_gt_count_cont,
                                        final_het_ab_cont=mt.index_rows(mt_orig.row_key).final_het_ab_cont)
        return mt_orig
    else:
        return mt


def filter_variants(mt, args, mode):
    """
    Low pass filter variants function, to run before samples QC to get rid of the 'worst' variants. Filters out GTs
    with low depth and bad quality, variants on call rate less than 80% (without correction for sex), snp and indel
    quality by depth, and variants whose fraction of het genotypes that are out of allelic balance is greater than 70%.

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
            hl.is_defined(mt.row.was_split)
        except Exception as e:
            logging.info("Split multi-allelics before running!")
            logging.info(e)
            return

        p_hwe = args.low_pass_p_hwe
        call_rate = args.low_pass_min_call_rate
        annotation_name = "variant_qc_thresholds_lowpass"
        skip_sex_handling = True

    elif mode == "final":
        # Check that samples QC has been run
        try:
            hl.is_defined(mt.failing_qc_samples)
            hl.is_defined(mt.pop_outlier_samples)
        except Exception as e:
            logging.info('failing_qc_samples and pop_outlier_samples not defined! Run samples QC before invoking this.')
            logging.info(e)

        p_hwe = args.final_p_hwe
        call_rate = args.final_min_call_rate
        annotation_name = "variant_qc_thresholds_final"
        skip_sex_handling = False

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

    annotation_dict = {'pval_hwe': str(p_hwe), 'hwe_excluded_in': str(hwe_tag),
                       'call_rate': str(call_rate),
                       'skipped_sex_handling': skip_sex_handling, 'snp_qd': str(args.snp_qd),
                       'indel_qd': str(args.indel_qd),
                       'het_max_ref_reads_thresh': str(args.max_het_ref_reads),
                       'het_min_ref_reads_thresh': str(args.min_het_ref_reads)}

    mt = mt.annotate_globals(**{annotation_name: annotation_dict})

    ###########################################################
    # Do genotype QC filtering and allelic balance annotation #
    ###########################################################
    if mode == "low_pass":
        # FIRST filter on genotype depth and quality, THEN calculate allelic balance across GTs for each variant
        mt = filter_genotypes_depth_quality(mt, args)

        # Define het allelic balance value (percentage of GT calls for that variant *in balance*)
        mt = annotate_ab(mt, args, prefix='lowpass_', filter_bad_samples=False)
    
    else:
        # Define het allelic balance value (percentage of GT calls for that variant *in balance*)
        mt = annotate_ab(mt, args, prefix='final_', filter_bad_samples=True)

        mt = mt.checkpoint(args.output_stem + '_filtervartemp1_deleteme.mt', overwrite=True)
        
        # Finally filter out those out of balance genotypes
        mt = filter_genotypes_ab(mt, args)

        #############################################################################
        # Filter out samples not passing samples QC, copy original mt for later use #
        #############################################################################
        mt_orig = mt
        mt = mt.filter_cols((mt.non_finns_to_remove == True) | (mt.fail_analytical > 0), keep=False)
        logging.info("Excluding samples failing qc. Number of samples on which variant QC is calculated: " +
                     str(mt.count_cols()))

        mt = mt.checkpoint(args.output_stem + '_filtervartemp2_deleteme.mt', overwrite=True)

        # Run variant QC after filtering out 'bad' samples, calculate call rate in hets
        logging.info('Running Hail variant QC function.')
        mt = hl.variant_qc(mt, name='final_no_failing_samples_varqc')
        mt = mt.annotate_rows(final_het_callrate=hl.agg.count_where(hl.is_defined(mt.GT) & mt.GT.is_het()) /
                                                 mt.final_het_gt_count)

        # Sanity check: is the het callrate ever higher than het allelic balance?
        logging.info('Any variants where het callrate is higher than het allelic balance?')
        logging.info(mt.aggregate_rows(hl.agg.count_where(mt.final_het_callrate > mt.final_het_ab_all)))

    ######################################################
    # Set HWE calculation in either all or only controls #
    ######################################################
    if args.pheno_col is not None:
        mt = mt.annotate_rows(
            hwe_ctrls_only=hl.agg.filter(((mt[args.pheno_col] == False) | hl.is_missing(mt.col[args.pheno_col])),
                                         hl.agg.hardy_weinberg_test(mt.GT)))
        logging.info('Annotated with hwe_ctrls_only successfully')
    else:
        logging.info('Phenotype column name not provided, calculating HWE for variants using all samples.')

    mt = mt.checkpoint(args.output_stem + '_filtervartemp3_deleteme.mt', overwrite=True)

    ###################
    # Filter variants #
    ###################

    # Initial count of variants
    varcount1 = mt.count_rows()
    logging.info('Finished initial variant count.')

    # Remove those with GSQR filters not empty (set length == 0)
    mt = mt.filter_rows(hl.len(mt.filters) != 0, keep=False)
    varcount2 = mt.count_rows()
    filter_val = varcount1 - varcount2
    logging.info('Finished variant count 2/')

    # Filter on HWE pvalue (controls only, or all)
    if args.pheno_col is not None:
        hwe_statement = mt.row.hwe_ctrls_only.p_value < args.low_pass_p_hwe
    else: hwe_statement = mt.row.lowpass_variant_qc.p_value_hwe < args.low_pass_p_hwe
    mt = mt.filter_rows(hwe_statement, keep=False)
    varcount3 = mt.count_rows()
    phwe_val = varcount2 - varcount3

    # remove those with less than min call rate
    mt = mt.filter_rows(mt.lowpass_variant_qc.call_rate < min_call_rate, keep=False)
    varcount4 = mt.count_rows()
    callrate_val = varcount3 - varcount4

    # If snp/indel, remove those with QD less than snp/indel QD threshold
    snps = hl.is_snp(mt.row.alleles[0], mt.row.alleles[1]) & (mt.row.info.QD < snp_qd)
    indels = (~hl.is_snp(mt.row.alleles[0], mt.row.alleles[1])) & (mt.row.info.QD < indel_qd)

    mt = mt.filter_rows((snps | indels), keep=False)
    varcount5 = mt.count_rows()
    snp_or_indel_val = varcount4 - varcount5

    varcount6 = mt.count_rows()
    # remove those with allele count = 0 (no variants exist in the population)
    mt = mt.filter_rows(hl.sum(mt.row.lowpass_variant_qc.AC) == hl.int(0), keep=False)
    ac_val = varcount5 - varcount6

    # Annotate exclusion counts to globals
    mt = mt.annotate_globals(variant_qc_exclusions_lowpass=
                                       {'mt.row.filters': filter_val, 'p_hwe': phwe_val, 'callrate': callrate_val,
                                        'sum_ac_0': ac_val, 'snp_or_indel_qd': snp_or_indel_val})

    # Report exclusion numbers to logs

    logging.info("\nThrowing out variants with length mt.row.filters != 0 (" + str(filter_val) + ")" +
                 ",\np value HWE < " + str(p_hwe) + " (" + str(phwe_val) + ")" +
                 " in " + hwe_tag +
                 "\ncall rate < " + str(min_call_rate) + " (" + str(callrate_val) + ")" +
                 ", \nAC = 0 (" + str(ac_val) + ")" +
                 ", \nsnp QD >= " + str(snp_qd) +
                 " or indel QD >= " + str(indel_qd) + " (" + str(snp_or_indel_val) + ")" )

    logging.info("\nAfter variant filter counts:" + str(mt.count_rows()))

    return mt