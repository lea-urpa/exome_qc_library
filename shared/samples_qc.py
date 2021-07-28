"""
Functions for use in samples QC for exome sequencing data with Hail.

Author: Lea Urpa, August 2020
"""
import logging
import time
import hail as hl
import utils
from bokeh.io import output_file, save
import networkx as nx


def filter_failing(mt, checkpoint_name, prefix="", pheno_col=None, entries=True, variants=True, samples=True,
                   unfilter_entries=False, pheno_qc=False, min_dp=10, min_gq=20, max_het_ref_reads=0.8,
                   min_het_ref_reads=0.2, min_hom_ref_ref_reads=0.9, max_hom_alt_ref_reads=0.1):
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
    if not prefix.endswith("_"):
        prefix = prefix + "_"

    if checkpoint_name.endswith(".mt"):
        checkpoint_name = checkpoint_name.replace(".mt", "")
    elif checkpoint_name.endswith(".mt/"):
        checkpoint_name = checkpoint_name.replace(".mt/", "")

    start_count = mt.count()
    tag = []

    ##################
    # Filter entries #
    ##################
    if entries:
        dp_cond = hl.is_defined(mt.DP) & (mt.DP > min_dp)
        gq_cond = hl.is_defined(mt.GQ) & (mt.GQ > min_gq)

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

        mt = mt.filter_entries(hl.is_defined(mt.GT) & dp_cond & gq_cond & (het_ab_cond | hom_ab_cond | homalt_ab_cond))
        mt = mt.checkpoint(checkpoint_name + "_GT_filtered.mt/", overwrite=True)
        tag.append("entries")

    ###################
    # Filter variants #
    ###################
    if variants:
        mt = mt.filter_rows((hl.len(mt[prefix + 'failing_variant_qc']) == 0) &
                            hl.is_defined(mt[prefix + "failing_variant_qc"]), keep=True)
        tag.append("variants")

        if (pheno_col is not None) and (pheno_qc is True):
            mt = mt.filter_rows((hl.len(mt.failing_pheno_varqc) == 0) & hl.is_defined(mt.failing_pheno_varqc),
                                keep=True)
            tag.append("variants by phenotype")

        mt = mt.checkpoint(checkpoint_name + "_variants_filtered.mt/", overwrite=True)

    ##################
    # Filter samples #
    ##################
    if samples:
        mt = mt.filter_cols((hl.len(mt.failing_samples_qc) == 0) & hl.is_defined(mt.failing_samples_qc) &
                            (mt.pop_outlier_sample == False) & hl.is_defined(mt.pop_outlier_sample), keep=True)

        mt = mt.checkpoint(checkpoint_name + "_samples_filtered.mt/", overwrite=True)
        tag.append("samples")

    final_count = mt.count()

    logging.info(f"Removing failing {', '.join(tag)}.")
    if entries:
        logging.info("(including entries failing allelic balance). "
                     "Not filtering genotypes missing DP, GQ, AD, or PL measures.")

    logging.info(f"Matrix table count before filtering: {start_count}. After filtering: {final_count}")

    if (unfilter_entries is True) and (entries is True):  # Unfilter entries if needed for pc_relate
        logging.info('Unfiltering entries (setting them to missing)')
        mt = mt.unfilter_entries()

    return mt


def find_pop_outliers(mt, checkpoint_name, pop_sd_threshold=4, plots=True, max_iter=8, reference_genome="GRCh38",
                      pca_plot_annotations=None):
    """
    Takes an LD pruned matrix table and determines which individuals are not clustering with Finns,
    with the principal components standard-deviation method, a la the FinnGen genotype team.
    Basically, calculates first two principal components, excludes outliers > 4 SD on those two
    components, recalculates, and repeats until there are no outliers.

    :param mt: matrix table from which to calculate principal components and find population outliers
    :param plots: plot the principal components for each round?
    :param max_iter: maximum number of iterations on running the PC plots, mostly for testing purposes.
    :return: returns the unfiltered, annotated matrix table
    """
    outliers = 1
    round_num = 1
    datestr = time.strftime("%Y.%m.%d")

    ##################
    # Filter dataset #
    ##################

    # Filter to autosomes
    if reference_genome is "GRCh38":
        autosomes = ["chr" + str(i) for i in range(1, 23)]
    else:
        autosomes = [str(i) for i in range(1, 23)]

    mt_autosomes = mt.filter_rows(hl.literal(autosomes).contains(mt.locus.contig))
    var_count = mt_autosomes.count_rows()
    logging.info(f"Variant count after filtering to autosomes: {var_count}")

    # Remove related individuals
    mt = mt.filter_cols(mt.related_to_remove == False)
    sample_count = mt.count_cols()
    logging.info(f'Sample count after removing related individuals: {sample_count}')

    # Write checkpoint
    mt_fn = checkpoint_name.rstrip("/").replace(".mt", "") + "_autosomes_norelatives"
    mt = mt.checkpoint(mt_fn)

    #########################################################################
    # Calculate PCAs, detect outliers, and repeat until outliers count == 0 #
    #########################################################################
    while outliers > 0:
        # Check if round is too big
        if round_num > max_iter:
            logging.warning('Warning- number of iterations in PC plots too many. Check the PC plots and see '
                            'what is wrong! Cutting rest of iterations, moving on in pipeline.')
            outliers = 0
            continue

        # Calculate pcs
        mt_count = mt_autosomes.count()
        logging.info(f"Count of samples and variants for matrix table PCA is calculated on: {mt_count}")
        eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt_autosomes.GT, k=2)

        # Parse column annotations for PCA plots
        if pca_plot_annotations is not None:
            pca_annotations = pca_plot_annotations.strip().split(",")

        # Do PCA plots, if specified
        if plots:
            coldata = mt.cols()
            if pca_plot_annotations is not None:
                for annotation in pca_annotations:
                    scores = scores.annotate(**{annotation: coldata[scores.s][annotation]})

            if pca_plot_annotations is not None:
                for annotation in pca_annotations:
                    output_file(f"{datestr}_find_population_outliers_pcsplots_round{round_num}_{annotation}.html")
                    p = hl.plot.scatter(scores.scores[0], scores.scores[1], label=scores[annotation])
                    save(p)
            else:
                output_file(f"{datestr}_find_population_outliers_pcsplots_round{round_num}.html")
                p = hl.plot.scatter(scores.scores[0], scores.scores[1])
                save(p)

        # Calculate upper and lower limits for each principal component
        pc1_stats = scores.aggregate(hl.agg.stats(scores.scores[0]))
        cutoff1_upper = pc1_stats.mean + (pop_sd_threshold * pc1_stats.stdev)
        cutoff1_lower = pc1_stats.mean - (pop_sd_threshold * pc1_stats.stdev)

        pc2_stats = scores.aggregate(hl.agg.stats(scores.scores[1]))
        cutoff2_upper = pc2_stats.mean + (pop_sd_threshold * pc2_stats.stdev)
        cutoff2_lower = pc2_stats.mean - (pop_sd_threshold * pc2_stats.stdev)

        # Make table of those failing cutoffs
        pc1_cut = (scores.scores[0] > cutoff1_upper) | (scores.scores[0] < cutoff1_lower)
        pc2_cut = (scores.scores[1] > cutoff2_upper) | (scores.scores[1] < cutoff2_lower)
        outlier_table = scores.filter(pc1_cut | pc2_cut, keep=True)

        # Get iterator for number of failing individuals
        outliers = outlier_table.count()
        logging.info(f"Number of individuals outside mean +/- {pop_sd_threshold} standard deviations for PCS 1 or "
                     f"2 in  round {round_num}: {outliers}")

        # Add round info to table and join to bigger table
        outlier_table = outlier_table.annotate(round=round_num)
        if round_num == 1:
            all_outliers = outlier_table
        else:
            all_outliers = all_outliers.union(outlier_table)

        # Filter outlier samples from main table to run PCA calculation again
        mt_autosomes = mt_autosomes.anti_join_cols(outlier_table)

        round_num += 1

    # Write table of all outliers
    pop_outliers_fn = checkpoint_name.rstrip("/").replace(".mt", "") + "_population_outliers.txt"
    all_outliers.export(pop_outliers_fn)
    pop_outliers = all_outliers.s.take(all_outliers.count())

    # Report total number of population outliers
    logging.info(f"Total number of population outliers: {len(pop_outliers)}")

    # Return mt and samples list
    return pop_outliers


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

    # Pull data to cols and checkpoint
    mt_cols = mt.cols()
    mt_cols = mt_cols.checkpoint("samples_qc_cols_tmp.ht", overwrite=True)

    # Instantiate empty array for failing samples QC tags
    mt_cols = mt_cols.annotate(failing_samples_qc=hl.empty_array(hl.tstr))

    ############################################################
    # Find samples failing on chimeras or contamination values #
    ############################################################
    mt_cols = mt_cols.annotate(failing_samples_qc=hl.cond(
        (mt_cols[args.chimeras_col] > args.chimeras_max) & hl.is_defined(mt_cols[args.chimeras_col]),
        mt_cols.failing_samples_qc.append("failing_chimeras"),
        mt_cols.failing_samples_qc))

    mt_cols = mt_cols.annotate(failing_samples_qc=hl.cond(
        (mt_cols[args.contamination_col] > args.contamination_max) & hl.is_defined(mt_cols[args.contamination_col]),
        mt_cols.failing_samples_qc.append("failing_contamination"),
        mt_cols.failing_samples_qc))

    failing_chim = mt_cols.aggregate(hl.agg.count_where(mt_cols.failing_samples_qc.contains("failing_chimeras")))
    miss_chim = mt_cols.aggregate(hl.agg.count_where(~(hl.is_defined(mt_cols[args.chimeras_col]))))
    failing_contam = mt_cols.aggregate(hl.agg.count_where(mt_cols.failing_samples_qc.contains("failing_contamination")))
    miss_contam = mt_cols.aggregate(hl.agg.count_where(~(hl.is_defined(mt_cols[args.contamination_col]))))

    logging.info(f"Number of samples failing on chimeras % > {args.chimeras_max}: {failing_chim}")
    logging.info(f"Number of samples missing chimeras %: {miss_chim}")
    logging.info(f"Number of samples failing on contamination % > {args.contamination_max}: {failing_contam}")
    logging.info(f"Number of samples missing contamination %: {miss_contam}")

    chim_stats = mt_cols.aggregate(hl.agg.stats(mt_cols[args.chimeras_col]))
    cont_stats = mt_cols.aggregate(hl.agg.stats(mt_cols[args.contamination_col]))
    logging.info(f"Chimeras statistics: {chim_stats}")
    logging.info(f"Contamination statistics: {cont_stats}")

    ###############################################
    # Find samples failing on sex-aware call rate #
    ###############################################
    if args.sample_call_rate is not None:
        mt_cols = mt_cols.annotate(failing_samples_qc=hl.cond(
            (mt_cols.sexaware_sample_call_rate < args.sample_call_rate) & hl.is_defined(mt_cols.sexaware_sample_call_rate),
            mt_cols.failing_samples_qc.append("failing_sexaware_sample_call_rate"),
            mt_cols.failing_samples_qc))

        mt_cols = mt_cols.annotate(failing_samples_qc=hl.cond(
            ~(hl.is_defined(mt_cols.sexaware_sample_call_rate)),
            mt_cols.failing_samples_qc.append("missing_sexaware_sample_call_rate"),
            mt_cols.failing_samples_qc))

        failing_cr = mt_cols.aggregate(hl.agg.count_where(mt_cols.failing_samples_qc.contains("failing_sexaware_sample_call_rate")))
        missing_cr = mt_cols.aggregate(
            hl.agg.count_where(mt_cols.failing_samples_qc.contains("missing_sexaware_sample_call_rate")))

        logging.info(f"Number of samples failing on sex-aware call rate > {args.sample_call_rate}: {failing_cr}")
        logging.info(f"Number of samples missing sex-aware call rate : {missing_cr}")

        cr_stats = mt_cols.aggregate(hl.agg.stats(mt_cols.sexaware_sample_call_rate))

        logging.info(f"Sex-aware call rate statistics: {cr_stats}")

    ######################################################################################
    # Find samples failing per-cohort on titv, het_homvar ratio, indel, and # singletons #
    ######################################################################################
    if args.batch_col_name is not None:
        batch_none = mt_cols.aggregate(hl.agg.count_where(~(hl.is_defined(mt_cols[args.batch_col_name]))))
        mt_cols = mt_cols.annotate(**{args.batch_col_name: hl.or_else(mt_cols[args.batch_col_name], "no_batch_info")})

        if batch_none > 0:
            logging.info(f"Warning- {batch_none} samples have batch undefined. These samples will be grouped in one"
                         f"batch for sample QC (named no_batch_info).")
            mt_cols.filter_cols(mt_cols[args.batch_col_name] == "no_batch_info").s.show(batch_none + 1)

        batch_set = mt_cols.aggregate(hl.agg.collect_as_set(mt_cols[args.batch_col_name]))
    else:
        args.batch_col_name = "mock_batch_col"
        mt_cols = mt_cols.annotate(mock_batch_col="all")
        batch_set = ["all"]

    # Convert batch strings to numeric values, create label for plotting
    batch_set_numeric = list(range(len(batch_set)))
    batch_key = list(zip(batch_set, batch_set_numeric))

    mt_cols = mt_cols.annotate(plot_batch=0)
    for batch in batch_key:
        mt_cols = mt_cols.annotate(plot_batch=hl.cond(mt_cols[args.batch_col_name] == batch[0],
                                                      batch[1], mt_cols.plot_batch))
        mt_cols = mt_cols.annotate(plot_batch_jitter=mt_cols.plot_batch + hl.rand_unif(-0.3, 0.3))

    batch_thresholds = {}
    batch_statistics = {}
    for measure in ['r_ti_tv', 'r_het_hom_var', 'r_insertion_deletion', 'n_singleton']:
        logging.info(f"Performing sample QC for measure {measure}")

        # Instantiate/reset box plot label
        mt_cols = mt_cols.annotate(boxplot_label=mt_cols[args.batch_col_name])

        batch_thresholds[measure] = {}
        batch_statistics[measure] = {}

        mt_cols = mt_cols.annotate(failing_samples_qc=hl.cond(
            ~(hl.is_defined(mt_cols.sample_qc[measure])),
            mt_cols.failing_samples_qc.append(f"missing_{measure}"),
            mt_cols.failing_samples_qc))

        for batch in batch_set:
            # See if values exist at all for all values
            defined_values = mt_cols.aggregate(hl.agg.count_where(hl.is_defined(mt_cols.sample_qc[measure])))

            if defined_values > 0:
                # Get mean and standard deviation for each measure, for each batch's samples
                stats = mt_cols.aggregate(hl.agg.filter(mt_cols[args.batch_col_name] == batch,
                                                        hl.agg.stats(mt_cols.sample_qc[measure])))

                # Get cutoffs for each measure
                cutoff_upper = stats.mean + (args.sampleqc_sd_threshold * stats.stdev)
                cutoff_lower = stats.mean - (args.sampleqc_sd_threshold * stats.stdev)

                if measure == "n_singleton":
                    logging.info(f"Max number of singletons for batch {batch}: {stats.max}")

                mt_cols = mt_cols.annotate(failing_samples_qc=hl.cond(
                    ((mt_cols.sample_qc[measure] > cutoff_upper) | (mt_cols.sample_qc[measure] < cutoff_lower))
                    & hl.is_defined(mt_cols.sample_qc[measure]) & (mt_cols[args.batch_col_name] == batch),
                    mt_cols.failing_samples_qc.append(f"failing_{measure}"),
                    mt_cols.failing_samples_qc))

                mt_cols = mt_cols.annotate(boxplot_label=hl.cond(
                    ((mt_cols.sample_qc[measure] > cutoff_upper) | (mt_cols.sample_qc[measure] < cutoff_lower)) &
                    hl.is_defined(mt_cols.sample_qc[measure]) & (mt_cols[args.batch_col_name] == batch),
                    "outlier", mt_cols.boxplot_label))

                # Collect thresholds and statistics for each batch
                batch_thresholds[measure][batch] = {'min_thresh': cutoff_lower, 'max_thresh': cutoff_upper}
                batch_statistics[measure][batch] = stats

            else:
                logging.error(f"Error- no defined values for measure {measure}. NAs can be introduced by division by "
                              f"zero. Samples not filtered on {measure}!")

        # Create plot for measure for each batch
        output_file(f"{datestr}_samples_qc_plots_{measure}.html")
        p = hl.plot.scatter(mt_cols.plot_batch_jitter, mt_cols.sample_qc[measure],
                            label=mt_cols.boxplot_label,
                            title=f"{measure} values split by batch.")
        save(p)


    ##########################
    # Report failing samples #
    ##########################
    for measure in ['r_ti_tv', 'r_het_hom_var', 'r_insertion_deletion', 'n_singleton']:
        failing_count = mt_cols.aggregate(hl.agg.count_where(mt_cols.failing_samples_qc.contains(f"failing_{measure}")))
        missing_count = mt_cols.aggregate(hl.agg.count_where(mt_cols.failing_samples_qc.contains(f"missing_{measure}")))
        logging.info(f"Number of samples failing on {measure}: {failing_count}")
        logging.info(f"Number of samples missing {measure}: {missing_count}")

    failing_any = mt_cols.aggregate(hl.agg.count_where(hl.len(mt_cols.failing_samples_qc) != 0))
    logging.info(f"Number of samples failing samples QC on any measure: {failing_any}")

    if args.pheno_col is not None:
        cases_failing = mt_cols.aggregate(hl.agg.filter(mt_cols[args.pheno_col] == True,
                                                           hl.agg.count_where(hl.len(mt_cols.failing_samples_qc) != 0)))
        controls_failing = mt_cols.aggregate(hl.agg.filter(mt_cols[args.pheno_col] == False,
                                                           hl.agg.count_where(hl.len(mt_cols.failing_samples_qc) != 0)))
        logging.info(f"Cases failing QC: {cases_failing}")
        logging.info(f"Controls failing QC: {controls_failing}")

    #######################################################################################################
    # Annotate original (unfiltered) matrix table with failing samples QC information + sample QC measure #
    #######################################################################################################
    mt_to_annotate = mt_to_annotate.annotate_cols(sample_qc=mt_cols[mt_to_annotate.s].sample_qc)
    mt_to_annotate = mt_to_annotate.annotate_cols(failing_samples_qc=mt_cols[mt_to_annotate.s].failing_samples_qc)

    mt_to_annotate = mt_to_annotate.annotate_globals(samples_qc_stats_batches=batch_statistics)
    mt_to_annotate = mt_to_annotate.annotate_globals(samples_qc_stats_chim_cont={'chimeras': chim_stats,
                                                                                 'contamination': cont_stats})
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

    output_file(f"{datestr}_imputed_sex_fstat_hist.html")
    p = hl.plot.histogram(fstat_hist, legend='F stat', title='F stat histogram')
    save(p)

    if mt_to_annotate is not None:
        mt_to_annotate = mt_to_annotate.annotate_cols(is_female_imputed=imputed_sex[mt_to_annotate.s].is_female,
                                                      f_stat=imputed_sex[mt_to_annotate.s].f_stat)
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


def filter_graph_cases(graph, cases):
    return [elem for elem in graph.nodes() if elem in cases]


def filter_node_cases(nodes, cases):
    return [elem for elem in nodes if elem in cases]


def sanity_check(g, nodes):
    """
    Given a list of nodes it makes sure that the algorithms are working properly.
    That is, that the subgraph induced by the remaining nodes does not contain edges.
    :param g: network x graph
    :param nodes: notes of a networkx graph
    :return:
    """
    assert g.subgraph(nodes).number_of_edges() == 0


def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c).copy()


def nx_algorithm(g, cases):
    """
    nx native based method for filtering cases. In each subgraph it maximizes the independent cases (read: disease
    affected) first and then proceeds with the rest of the subgraph (by including these cases as required nodes in the
    subsequent maximally independent graph)

    :param g: nx graph
    :param cases: list of cases
    :return: list of related nodes that are discarded
    """
    #####################################################
    # Report the number of highly connected individuals #
    #####################################################
    degrees = dict(g.degree())
    high_degree_count = 0
    for degree in degrees.values():
        if degree > 100:
            high_degree_count += 1

    if high_degree_count > 0:
        print(f"Warning! {high_degree_count} samples are connected to > 100 other samples. This is likely from poorly "
              f"calculated kinships, increase number of variants used to calculate kinships and try again.")

    ##########################################
    # Loop through subgraphs in larger graph #
    ##########################################
    # Instantiate list of unrelated notes
    unrelated_nodes = []
    num_graphs = 0
    for subgraph in connected_component_subgraphs(g):
        num_graphs += 1
        if cases is not None:
            subcases = filter_graph_cases(subgraph, cases)  # list of local cases
            if len(subcases) > 0:
                # graph induced by local cases
                case_graph = subgraph.subgraph(cases)
                # get maximal set of cases in case graph
                unrelated_subcases = nx.maximal_independent_set(case_graph)
                # get maximal set of nodes in subgraph, giving unrelated cases to keep
                unrelated_nodes += nx.maximal_independent_set(subgraph, unrelated_subcases)
            else:
                unrelated_nodes += nx.maximal_independent_set(subgraph)
        else:
            unrelated_nodes += nx.maximal_independent_set(subgraph)

    logging.info(f"Number of groups of related individuals: {num_graphs}")

    ####################
    # result summaries #
    ####################
    related_nodes = list(set(g.nodes()) - set(unrelated_nodes))
    unrelated_cases = filter_node_cases(unrelated_nodes, cases)

    #####################################################################################
    # sanity_check: makes sure that the unrelated cases/nodes are in fact not connected #
    #####################################################################################
    sanity_check(g, unrelated_nodes)
    sanity_check(g, unrelated_cases)

    return related_nodes, len(unrelated_cases), len(unrelated_nodes)


def king_relatedness(mt, checkpoint_name, kinship_threshold=0.0883, pheno_col=None, force=False,
                     cluster_name=None, num_secondary_workers=None, region=None, reference_genome="GRCh38"):
    """

    :param mt: matrix table to calculate kinship values from
    :param checkpoint_name: checkpoint name (including bucket, if applicable)
    :param kinship_threshold: kinship threshold to declare related pairs
    :param pheno_col: column in matrix table (boolean) indicating case status
    :param force: force re-analysis of checkpointed steps?
    :param cluster_name: name of cluster to update with secondary workers during kinship calculation
    :param num_secondary_workers: number of secondary workers to add
    :param region: region of cluster
    :param reference_genome: reference genome to use (for filtering to autosomes only)
    :return:
    """

    kinship_fn = checkpoint_name.rstrip("/").replace(".mt", "") + "_kinship.mt/"
    relatives_fn = checkpoint_name.rstrip("/").replace(".mt", "") + "_related_pairs.txt"

    # Filter to autosomes
    if reference_genome is "GRCh38":
        autosomes = ["chr" + str(i) for i in range(1, 23)]
    else:
        autosomes = [str(i) for i in range(1,23)]
    mt_autosomes = mt.filter_rows(hl.literal(autosomes).contains(mt.locus.contig))

    var_count = mt_autosomes.count_rows()
    if var_count < 10000:
        logging.warning("Warning! Number of variants to calculate kinship is less than 10 000, it is likely that "
                        "kinship calculations will be incorrect (more relatedness than reality). Number of variants "
                        f"left after removing sex chromosomes: {var_count}.")

    # Calculate kinship
    if (not utils.check_exists(kinship_fn)) or force:
        if (cluster_name is not None) and (num_secondary_workers is not None) and (region is not None):
            utils.add_secondary(cluster_name, num_secondary_workers, region)
        kinship = hl.king(mt_autosomes.GT)
        logging.info(f"Writing kinship matrix table to file: {kinship_fn}")
        kinship = kinship.checkpoint(kinship_fn, overwrite=True)
        if (cluster_name is not None) and (num_secondary_workers is not None) and (region is not None):
            utils.remove_secondary(cluster_name, region)
    else:
        logging.info(f"Detected kinship file exists {kinship_fn}, loading that.")
        kinship = hl.read_matrix_table(kinship_fn)

    # Get just pairs above threshold, convert to pandas df
    relatives = kinship.filter_entries(kinship.phi > kinship_threshold)

    rel_tab = relatives.entries()
    rel_tab = rel_tab.filter(rel_tab.s != rel_tab.s_1)  # remove diagonal
    logging.info(f"Writing kinship pairs above threshold table to file: {relatives_fn}")
    rel_tab.export(relatives_fn)

    # Select just edges and convert to pandas df
    rel_tab = rel_tab.key_by().select("s", "s_1")
    rel_df = rel_tab.to_pandas()

    # Create graph from pandas df
    related_ind_g = nx.from_pandas_edgelist(rel_df, "s", "s_1")

    if pheno_col is not None:
        sample_info = mt.cols()
        case_info = sample_info.filter(sample_info.is_case == True)
        case_ids = case_info.s.take(case_info.count())
    else:
        case_ids = None

    logging.info('Calculating maximal independent set.')
    related_to_remove, num_ind_cases, num_ind_nodes = nx_algorithm(related_ind_g, case_ids)
    logging.info('# of unrelated cases given king input: ' + str(num_ind_cases))

    # Annotate matrix table with maximal unrelated individuals list
    mt = mt.annotate_cols(related_to_remove=hl.if_else(hl.literal(related_to_remove).contains(mt.s), True, False))
    return mt, related_to_remove
