"""
Functions for use in samples QC for exome sequencing data with Hail.

Author: Lea Urpa, August 2020
"""
import logging
import time
import hail as hl
from bokeh.io import output_file, save


def filter_failing(mt, args, mode, entries=True, variants=True, samples=True, unfilter_entries=False, pheno_qc=False):
    """
    Filters failing samples, variants, and entries from a given matrix table
    :param mt: matrix table to filter
    :param args:
    :param mode: Filter on final or low_pass variant QC?
    :param entries: filter entries?
    :param variants: filter variants?
    :param samples: filter samples?
    :param unfilter_entries: Unfilter entries (set to missing) after filtering entries?
    :return:
    """
    start_count = mt.count()
    tag = []
    if entries:
        mt = mt.filter_entries((hl.len(mt[mode + '_failing_depth_quality']) == 0) &
                                hl.is_defined(mt[mode + '_failing_depth_quality']) &
                               (hl.len(mt[mode + '_failing_ab']) == 0) & hl.is_defined(mt[mode + "_failing_ab"]),
                               keep=True)
        tag.append("entries")
    if variants:
        mt = mt.filter_rows((hl.len(mt[mode + '_failing_variant_qc']) == 0) &
                            hl.is_defined(mt[mode + "_failing_variant_qc"]), keep=True)
        tag.append("variants")
        if (args.pheno_col is not None) and (pheno_qc is True):
            mt = mt.filter_rows((hl.len(mt.failing_pheno_varqc) == 0) & hl.is_defined(mt.failing_pheno_varqc),
                                keep=True)
            tag.append("variants by phenotype")

    if samples:
        mt = mt.filter_cols((hl.len(mt.failing_samples_qc) == 0) & hl.is_defined(mt.failing_samples_qc) &
                            (mt.pop_outlier_sample == False) & hl.is_defined(mt.pop_outlier_sample), keep=True)
        tag.append("samples")

    final_count = mt.count()

    logging.info(f"Removing failing {', '.join(tag)}.")
    if entries:
        logging.info("(including entries failing allelic balance)")

    logging.info(f"Matrix table count before filtering: {start_count}. After filtering: {final_count}")

    if unfilter_entries is True:  # Unfilter entries if needed for pc_relate
        logging.info('Unfiltering entries (setting them to missing)')
        mt = mt.unfilter_entries()

    logging.info(f"Writing temporary checkpoint for filtered mt.")
    mt = mt.checkpoint(f"{args.output_stem}_{mode}_removed_tmp_{args.tmp_counter}.mt", overwrite=True)
    args.tmp_counter = args.tmp_counter + 1

    return mt


def find_pop_outliers(mt_ldpruned, mt_to_annotate, args, plots=True, max_iter=8):
    """
    Takes an LD pruned matrix table and determines which individuals are not clustering with Finns,
    with the principal components standard-deviation method, a la the FinnGen genotype team.
    Basically, calculates first two principal components, excludes outliers > 4 SD on those two
    components, recalculates, and repeats until there are no outliers.

    :param mt_ldpruned: LD pruned matrix table from which to calculate principal components. Relatives should already
    be removed!
    :param mt_to_annotate:  Matrix table with full variants to annotate which individuals are outliers.
    :param plots: plot the principal components for each round?
    :param max_iter: maximum number of iterations on running the PC plots, mostly for testing purposes.
    :return: returns the unfiltered, annotated matrix table
    """
    # Instantiate counter
    outliers = 1
    i = 1
    datestr = time.strftime("%Y.%m.%d")

    # Overwrite pop_outlier_sample annotation in mt to annotate
    logging.info("Overwriting any previous information in column annotation pop_outlier_sample in mt to annotate.")
    mt_to_annotate = mt_to_annotate.annotate_cols(pop_outlier_sample=False)

    # Remove related individuals for population outlier analysis
    logging.info('Removing related individuals from LD pruned dataset for population outlier analysis '
                 '(keeping in main dataset).')
    mt_ldpruned = mt_ldpruned.filter_cols(mt_ldpruned.related_to_remove == True, keep=False)
    logging.info(f'Sample count after removing related individuals: {mt_ldpruned.count_cols()}')

    # Remove chromosome X
    if args.reference_genome == "GRCh38":
        chrom = "chrX"
    elif args.reference_genome == "GRCh37":
        chrom = "X"
    mt_ldpruned = mt_ldpruned.filter_rows(mt_ldpruned.locus.contig == chrom, keep=False)
    logging.info(f"Variant count after removing chromosome x: {mt_ldpruned.count_rows()}")

    #########################################################################
    # Calculate PCAs, detect outliers, and repeat until outliers count == 0 #
    #########################################################################
    while outliers > 0:
        # Check if i is too big
        if i > max_iter:
            logging.warning('Warning- number of iterations in PC plots too many. Check the PC plots and see '
                            'what is wrong! Cutting rest of iterations, moving on in pipeline.')
            outliers = 0
            continue

        # Calculate pcs
        logging.info(f"Count of samples and variants for matrix table PCA is calculated on: {mt_ldpruned.count()}")
        eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt_ldpruned.GT, k=2)

        # Parse column annotations for PCA plots
        if args.pca_plot_annotations is not None:
            pca_annotations = args.pca_plot_annotations.strip().split(",")

        # Do PCA plots, if specified
        if plots:
            coldata = mt_ldpruned.cols()
            if args.pca_plot_annotations is not None:
                for annotation in pca_annotations:
                    scores = scores.annotate(**{annotation: coldata[scores.s][annotation]})

            output_file(f"find_population_outliers_pcsplots_round{i}_{datestr}.html")
            if args.pca_plot_annotations is not None:
                for annotation in pca_annotations:
                    p = hl.plot.scatter(scores.scores[0], scores.scores[1], label=scores[annotation])
                    save(p)
            else:
                p = hl.plot.scatter(scores.scores[0], scores.scores[1])
                save(p)

        # Calculate upper and lower limits for each principal component
        pc1_stats = scores.aggregate(hl.agg.stats(scores.scores[0]))
        cutoff1_upper = pc1_stats.mean + (args.pop_sd_threshold * pc1_stats.stdev)
        cutoff1_lower = pc1_stats.mean - (args.pop_sd_threshold * pc1_stats.stdev)

        pc2_stats = scores.aggregate(hl.agg.stats(scores.scores[1]))
        cutoff2_upper = pc2_stats.mean + (args.pop_sd_threshold * pc2_stats.stdev)
        cutoff2_lower = pc2_stats.mean - (args.pop_sd_threshold * pc2_stats.stdev)

        # Make table of those failing cutoffs
        pc1_cut = (scores.scores[0] > cutoff1_upper) | (scores.scores[0] < cutoff1_lower)
        pc2_cut = (scores.scores[1] > cutoff2_upper) | (scores.scores[1] < cutoff2_lower)
        outlier_table = scores.filter(pc1_cut | pc2_cut, keep=True)

        # Get iterator for number of failing individuals
        outliers = outlier_table.count()
        logging.info(f"Number of individuals outside mean +/- {args.pop_sd_threshold} standard deviations for PCS 1 or "
                     f"2 in  round {i}: {outliers}")

        # Filter outlier samples from main table to run PCA calculation again
        mt_ldpruned = mt_ldpruned.anti_join_cols(outlier_table)

        # Annotate main matrix table with outlier samples (if sample exists in outlier_table, mark as outlier in mt)
        mt_to_annotate = mt_to_annotate.annotate_cols(pop_outlier_sample=
                                                      hl.cond(hl.is_defined(outlier_table[mt_to_annotate.s].scores),
                                                              True, mt_to_annotate.pop_outlier_sample))
        i = i + 1

    # Report total number of population outliers
    total_outliers = mt_to_annotate.aggregate_cols(hl.agg.count_where(mt_to_annotate.pop_outlier_sample == True))
    logging.info(f"Total number of population outliers: {total_outliers}")

    logging.info('Samples identified as population outliers:')
    mt_to_annotate.filter_cols(mt_to_annotate.pop_outlier_sample == True, keep=True).s.show()

    # Return mt and samples list
    return mt_to_annotate


def samples_qc(mt, mt_to_annotate, args):
    """
    Performs samples QC on a matrix table, removing samples on chimera and contamination %, as well as being +/- 4
    standard deviations from mean on TiTv, het/homvar, insertion/deletion ratios and n_singletons for a specific
    batch or cohort

    :param mt: matrix table, low-pass failing variants and genotypes filtered out
    :param mt_to_annotate: matrix table to annotate with failing samples information after calculating on filtered mt
    :param args:
    :return: returns annotated, unfiltered matrix table
    """
    datestr = time.strftime("%Y.%m.%d")

    # Run variant QC to get up to date variant QC metrics for samples QC
    mt = hl.sample_qc(mt)

    # Instantiate empty array for failing samples QC tags
    mt = mt.annotate_cols(failing_samples_qc=hl.empty_array(hl.tstr))

    ############################################################
    # Find samples failing on chimeras or contamination values #
    ############################################################
    mt = mt.annotate_cols(failing_samples_qc=hl.cond(
        (mt[args.chimeras_col] > args.chimeras_max) & hl.is_defined(mt[args.chimeras_col]),
        mt.failing_samples_qc.append("failing_chimeras"),
        mt.failing_samples_qc))

    mt = mt.annotate_cols(failing_samples_qc=hl.cond(
        ~(hl.is_defined(mt[args.chimeras_col])),
        mt.failing_samples_qc.append("missing_chimeras"),
        mt.failing_samples_qc))

    mt = mt.annotate_cols(failing_samples_qc=hl.cond(
        (mt[args.contamination_col] > args.contamination_max) & hl.is_defined(mt[args.contamination_col]),
        mt.failing_samples_qc.append("failing_contamination"),
        mt.failing_samples_qc))

    mt = mt.annotate_cols(failing_samples_qc=hl.cond(
        ~(hl.is_defined(mt[args.contamination_col])),
        mt.failing_samples_qc.append("missing_contamination"),
        mt.failing_samples_qc))

    failing_chim = mt.aggregate_cols(hl.agg.count_where(mt.failing_samples_qc.contains("failing_chimeras")))
    miss_chim = mt.aggregate_cols(hl.agg.count_where(mt.failing_samples_qc.contains("missing_chimeras")))
    failing_contam = mt.aggregate_cols(hl.agg.count_where(mt.failing_samples_qc.contains("failing_contamination")))
    miss_contam = mt.aggregate_cols(hl.agg.count_where(mt.failing_samples_qc.contains("missing_contamination")))

    logging.info(f"Number of samples failing on chimeras % > {args.chimeras_max}: {failing_chim}")
    logging.info(f"Number of samples missing chimeras %: {miss_chim}")
    logging.info(f"Number of samples failing on contamination % > {args.contamination_max}: {failing_contam}")
    logging.info(f"Number of samples missing contamination %: {miss_contam}")

    ###############################################
    # Find samples failing on sex-aware call rate #
    ###############################################
    mt = mt.annotate_cols(failing_samples_qc=hl.cond(
        (mt.sexaware_sample_call_rate < args.sample_call_rate) & hl.is_defined(mt.sexaware_sample_call_rate),
        mt.failing_samples_qc.append("failing_sexaware_sample_call_rate"),
        mt.failing_samples_qc))

    mt = mt.annotate_cols(failing_samples_qc=hl.cond(
        ~(hl.is_defined(mt.sexaware_sample_call_rate)),
        mt.failing_samples_qc.append("missing_sexaware_sample_call_rate"),
        mt.failing_samples_qc))

    failing_cr = mt.aggregate_cols(hl.agg.count_where(mt.failing_samples_qc.contains("failing_sexaware_sample_call_rate")))
    missing_cr = mt.aggregate_cols(
        hl.agg.count_where(mt.failing_samples_qc.contains("missing_sexaware_sample_call_rate")))

    logging.info(f"Number of samples failing on sex-aware call rate > {args.sample_call_rate}: {failing_cr}")
    logging.info(f"Number of samples missing sex-aware call rate : {missing_cr}")

    cr_stats = mt.aggregate_cols(hl.agg.stats(mt.sexaware_sample_call_rate))
    chim_stats = mt.aggregate_cols(hl.agg.stats(mt[args.chimeras_col]))
    cont_stats = mt.aggregate_cols(hl.agg.stats(mt[args.contamination_col]))

    logging.info(f"Sex-aware call rate statistics: {cr_stats}")
    logging.info(f"Chimeras statistics: {chim_stats}")
    logging.info(f"Contamination statistics: {cont_stats}")

    ######################################################################################
    # Find samples failing per-cohort on titv, het_homvar ratio, indel, and # singletons #
    ######################################################################################
    if args.batch_col_name is not None:
        batch_set = mt.aggregate_cols(hl.agg.collect_as_set(mt[args.batch_col_name]))
    else:
        args.batch_col_name = "mock_batch_col"
        mt = mt.annotate_cols(mock_batch_col="all")
        batch_set = ["all"]

    batch_thresholds = {}
    for batch in batch_set:
        batch_thresholds[batch] = {}
        for measure in ['r_ti_tv', 'r_het_hom_var', 'r_insertion_deletion', 'n_singleton']:
            defined_values = mt.aggregate_cols(hl.agg.count_where(hl.is_defined(mt.sample_qc[measure])))

            if defined_values > 0:
                # Get mean and standard deviation for each measure, for each batch's samples
                stats = mt.aggregate_cols(hl.agg.filter(mt[args.batch_col_name] == batch,
                                                        hl.agg.stats(mt.sample_qc[measure])))

                # Get cutoffs for each measure
                cutoff_upper = stats.mean + (args.sampleqc_sd_threshold * stats.stdev)
                cutoff_lower = stats.mean - (args.sampleqc_sd_threshold * stats.stdev)

                if measure == "n_singleton":
                    logging.info(f"Max number of singletons for batch {batch}: {stats.max}")

                mt = mt.annotate_cols(failing_samples_qc=hl.cond(
                    (mt.sample_qc[measure] > cutoff_upper) & hl.is_defined(mt.sample_qc[measure]),
                    mt.failing_samples_qc.append(f"failing_{measure}"),
                    mt.failing_samples_qc))

                mt = mt.annotate_cols(failing_samples_qc=hl.cond(
                    (mt.sample_qc[measure] < cutoff_lower) & hl.is_defined(mt.sample_qc[measure]),
                    mt.failing_samples_qc.append(f"failing_{measure}"),
                    mt.failing_samples_qc))

                mt = mt.annotate_cols(failing_samples_qc=hl.cond(
                    ~(hl.is_defined(mt.sample_qc[measure])),
                    mt.failing_samples_qc.append(f"missing_{measure}"),
                    mt.failing_samples_qc))

                # Collect thresholds for each batch
                batch_thresholds[batch][measure] = {'min_thresh': cutoff_lower, 'max_thresh': cutoff_upper}
            else:
                logging.error(f"Error- no defined values for measure {measure}. NAs can be introduced by division by "
                              f"zero. Samples not filtered on {measure}!")

    ##########################################
    # Create box plots for samples QC values #
    ##########################################
    # Convert batch strings to numeric values
    batch_set_numeric = list(range(len(batch_set)))
    batch_key = list(zip(batch_set, batch_set_numeric))

    mt_cols = mt.cols()
    mt_cols = mt_cols.annotate(plot_batch=0)

    for batch in batch_key:
        mt_cols = mt_cols.annotate(plot_batch=hl.cond(mt_cols[args.batch_col_name] == batch[0],
                                                      batch[1], mt_cols.plot_batch))

    output_file(f"samples_qc_measure_plots_{datestr}.html")
    for measure in ['r_ti_tv', 'r_het_hom_var', 'r_insertion_deletion', 'n_singleton']:
        defined_values = mt.aggregate_cols(hl.agg.count_where(hl.is_defined(mt.sample_qc[measure])))
        if defined_values > 0:
            p = hl.plot.scatter(mt_cols.plot_batch, mt_cols.sample_qc[measure], label=mt_cols[args.batch_col_name],
                                title=f"{measure} values split by batch.")
            save(p)

    ##########################
    # Report failing samples #
    ##########################
    for measure in ['r_ti_tv', 'r_het_hom_var', 'r_insertion_deletion', 'n_singleton']:
        failing_count = mt.aggregate_cols(hl.agg.count_where(mt.failing_samples_qc.contains(f"failing_{measure}")))
        missing_count = mt.aggregate_cols(hl.agg.count_where(mt.failing_samples_qc.contains(f"missing_{measure}")))
        logging.info(f"Number of samples failing on {measure}: {failing_count}")
        logging.info(f"Number of samples missing {measure}: {missing_count}")

    failing_any = mt.aggregate_cols(hl.agg.count_where(hl.len(mt.failing_samples_qc) != 0))
    logging.info(f"Number of samples failing samples QC on any measure: {failing_any}")

    if args.pheno_col is not None:
        cases_failing = mt.aggregate_cols(hl.agg.filter(mt[args.pheno_col] == True,
                                                           hl.agg.count_where(hl.len(mt.failing_samples_qc) != 0)))
        controls_failing = mt.aggregate_cols(hl.agg.filter(mt[args.pheno_col] == False,
                                                           hl.agg.count_where(hl.len(mt.failing_samples_qc) != 0)))
        logging.info(f"Cases failing QC: {cases_failing}")
        logging.info(f"Controls failing QC: {controls_failing}")

    #######################################################################################################
    # Annotate original (unfiltered) matrix table with failing samples QC information + sample QC measure #
    #######################################################################################################
    mt_to_annotate = mt_to_annotate.annotate_cols(sample_qc=mt_cols[mt_to_annotate.s].sample_qc)
    mt_to_annotate = mt_to_annotate.annotate_cols(failing_samples_qc=mt_cols[mt_to_annotate.s].failing_samples_qc)

    mt_to_annotate = mt_to_annotate.annotate_globals(samples_qc_thresholds=
                                                     {'chimeras_max': str(args.chimeras_max),
                                                      'contamination_max': str(args.contamination_max),
                                                      'deviation_multiplier_threshold': str(args.sampleqc_sd_threshold),
                                                      'batches': str(batch_set),
                                                      'batch_cohort_name': str(args.batch_col_name)})

    mt_to_annotate = mt_to_annotate.annotate_globals(samples_qc_batch_thresholds=batch_thresholds)

    return mt_to_annotate


def impute_sex_plot(mt, args, mt_to_annotate=None):
    """
    Impute sex of individuals and plot resultant f stat values
    :param mt: maf pruned matrix table to caculate f stat values
    :param mt_to_annotate: matrix table to add sex information to
    :return: returns either annotated matrix table and imputed sex Hail table, if mt_to_annotate is not None,
    or else just the imputed sex Hail table.
    """
    datestr = time.strftime("%Y.%m.%d")
    imputed_sex = hl.impute_sex(mt.GT, female_threshold=args.female_threshold, male_threshold=args.male_threshold)

    sex_count = imputed_sex.aggregate(hl.agg.counter(imputed_sex.is_female))

    logging.info(f'Imputed sex count: {sex_count}')

    fstat_stats = imputed_sex.aggregate(hl.agg.stats(imputed_sex.f_stat))
    fstat_hist = imputed_sex.aggregate(hl.agg.hist(imputed_sex.f_stat, fstat_stats.min, fstat_stats.max, 50))

    output_file(f"Imputed_sex_fstat_hist_{datestr}.html")
    p = hl.plot.histogram(fstat_hist, legend='F stat', title='F stat histogram')
    save(p)

    if mt_to_annotate is not None:
        mt_to_annotate = mt_to_annotate.annotate_cols(is_female_imputed=imputed_sex[mt_to_annotate.s].is_female)
        mt_to_annotate = mt_to_annotate.annotate_globals(sex_imputation_thresholds=
                                                         {'female_threshold': args.female_threshold,
                                                          'male_threshold': args.male_threshold})

        mt = mt.annotate_cols(is_female_imputed=imputed_sex[mt.s].is_female)
        mt = mt.annotate_globals(sex_imputation_thresholds={'female_threshold': args.female_threshold,
                                                            'male_threshold': args.male_threshold})
        args.sex_col = "is_female_imputed"
        args.male_tag = False
        args.female_tag = True

        return mt, imputed_sex, mt_to_annotate
    else:
        return mt, imputed_sex


def pc_project(mt, loadings_ht, loading_location="loadings", af_location="pca_af"):
    """
    Projects samples in `mt` on pre-computed PCs.
    :param MatrixTable mt: MT containing the samples to project into previously calculated PCs
    :param Table loadings_ht: HT containing the PCA loadings and allele frequencies used for the PCA
    :param str loading_location: Location of expression for loadings in `loadings_ht`
    :param str af_location: Location of expression for allele frequency in `loadings_ht`
    :return: Hail Table with scores calculated from loadings in column `scores`
    :rtype: Table

    From Konrad Karczewski
    """
    n_variants = loadings_ht.count()

    # Annotate matrix table with pca loadings and af from other dataset which pcs were calculated from
    mt = mt.annotate_rows(
        pca_loadings=loadings_ht[mt.row_key][loading_location],
        pca_af=loadings_ht[mt.row_key][af_location]
    )

    # Filter to rows where pca_loadings and af are defined, and af > 0 and < 1
    mt = mt.filter_rows(hl.is_defined(mt.pca_loadings) & hl.is_defined(mt.pca_af) &
                        (mt.pca_af > 0) & (mt.pca_af < 1))

    # Calculate genotype normalization constant
    # Basically, mean centers and normalizes the genotypes under the binomial distribution so that they can be
    # multiplied by the PC loadings to get the projected principal components
    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt.pca_af) / hl.sqrt(n_variants * 2 * mt.pca_af * (1 - mt.pca_af))

    mt = mt.annotate_cols(scores=hl.agg.array_sum(mt.pca_loadings * gt_norm))

    return mt.cols().select('scores')


def project_pcs_relateds(mt_ldpruned, mt, covar_pc_num):
    """
    Tales LD pruned matrix table, calculates PCs, and projects those PCs back to related individuals included in mt
    :param mt_ldpruned: matrix table with relatives removed, maf and ld pruned
    :param mt: matrix table with relatives included
    :param covar_pc_num: Number of principal components as covariates to calculate
    :return: returns matrix table with relatives, with PCs annotated
    """
    logging.info('Calculating principal components, annotating main dataset.')
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt_ldpruned.GT, k=covar_pc_num, compute_loadings=True)

    # Project PCs to related individuals
    # mt of related individuals only, not pop outliers or failing samples QC
    related_mt = mt.filter_cols((mt.related_to_remove == True) & (mt.pop_outlier_sample == False) &
                                (hl.len(mt.failing_samples_qc) == 0), keep=True)
    mt_ldpruned = mt_ldpruned.annotate_rows(pca_af=hl.agg.mean(mt_ldpruned.GT.n_alt_alleles()) / 2)
    mtrows = mt_ldpruned.rows()
    loadings = loadings.annotate(pca_af=mtrows[loadings.locus, loadings.alleles].pca_af)
    related_scores = pc_project(related_mt, loadings)

    # Add pcs as annotations to main table
    mt = mt.annotate_cols(**{'pc' + str(k+1): scores[mt.s].scores[k]
                             for k in range(covar_pc_num)})
    # Explanation: for k principal components in range 0 to covar_pc_num-1,
    # make pc k+1 (to start at pc1 instead of pc0) be the corresponding score (keyed by mt.s) from the table scores

    # Add pcs for related individuals
    mt = mt.annotate_cols(**{'pc' + str(k+1): hl.or_else(mt['pc'+str(k+1)], related_scores[mt.s].scores[k])
                             for k in range(covar_pc_num)})
    # Explanation: for k principal components in range from 0 to (covar_pc_num-1)
    # give either the existing pcX, or if missing give the corresponding score (keyed by mt.s)
    # from the table related_scores

    return mt
