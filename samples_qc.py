"""
Functions for use in samples QC for exome sequencing data with Hail.

Author: Lea Urpa, August 2020
"""
import logging
import time
import numpy as np
import hail as hl
from bokeh.io import output_file, save


def filter_failing(mt, args, varqc_name, entries=True, variants=True, samples=True, unfilter_entries=False):
    """

    :param mt:
    :param args:
    :param varqc_name:
    :param unfilter_entries: Unfilter entries following maf pruning? Necessary to run PCA later.
    :return:
    """
    start_count = mt.count()
    tag = []
    if entries:
        mt = mt.filter_entries((hl.len(mt.failing_depth_quality) == 0) & hl.len(mt.failing_ab) == 0, keep=True)
        tag.append("entries")
    if variants:
        mt = mt.filter_variants(hl.len(mt[varqc_name]) == 0, keep=True)
        tag.append("variants")
    if samples:
        mt = mt.filter_samples((hl.len(mt.failing_samples_qc) == 0) & mt.population_outlier == False, keep=True)
        tag.append("samples")

    final_count = mt.count()

    logging.info(f"Removing failing {', '.join(tag)}.")
    if entries:
        logging.info("Including entries failing allelic balance.")

    logging.info(f"Matrix table count before filtering: {start_count}. After filtering: {final_count}")

    if unfilter_entries is True:  # Unfilter entries if needed for pc_relate
        logging.info('Unfiltering entries (setting them to missing)')
        mt = mt.unfilter_entries()

    logging.info(f"Writing temporary checkpoint for filtered mt.")
    mt = mt.checkpoint(f"{args.outputstem}_{varqc_name}_removed_tmp.mt")

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
    mt_ldpruned = hl.filter_rows(mt_ldpruned.locus.contig == chrom, keep=False)
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
    datestr = time.strftime("%Y.%m.%d")

    # Run variant QC to get up to date variant QC metrics for samples QC
    mt = hl.sample_qc(mt)

    # Instantiate empty array for failing samples QC tags
    mt = mt.annotate_cols(failing_samples_qc=hl.empty_array(hl.tstr))

    ############################################################
    # Find samples failing on chimeras or contamination values #
    ############################################################
    mt = mt.annotate_cols(failing_samples_qc=hl.cond(hl.float(mt[args.chimeras_col]) > args.chimeras_max,
                                                     mt.failing_samples_qc.append("failing_chimeras"),
                                                     mt.failing_samples_qc))

    mt = mt.annotate_cols(failing_samples_qc=hl.cond(hl.float(mt[args.contamination_col]) > args.contamination_max,
                                                     mt.failing_samples_qc.append("failing_contamination"),
                                                     mt.failing_samples_qc))

    failing_chim = mt.aggregate_cols(mt.failing_samples_qc.contains("failing_chimeras"))
    failing_contam = mt.aggregate_cols(mt.failing_samples_qc.contains("failing_contamination"))

    logging.info(f"Number of samples failing on chimeras %: {failing_chim}")
    logging.info(f"Number of samples failing on contamination %: {failing_contam}")

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
            # Get mean and standard deviation for each measure, for each batch's samples
            stats = mt.aggregate_cols(hl.agg.filter(mt[args.batch_col_name] == batch,
                                                    hl.agg.stats(mt.sample_qc[measure])))

            # Get cutoffs for each measure
            cutoff_upper = stats.mean + (args.sampleqc_sd_threshold * stats.stdev)
            cutoff_lower = stats.mean - (args.sampleqc_sd_threshold * stats.stdev)

            if measure == "n_singleton":
                logging.info(f"Max number of singletons for batch {batch}: {stats.max}")

            mt = mt.annotate_cols(failing_samples_qc=hl.cond(mt.sample_qc[measure] > cutoff_upper,
                                                             mt.failing_samples_qc.append(f"failing_{measure}"),
                                                             mt.failing_samples_qc))

            mt = mt.annotate_cols(failing_samples_qc=hl.cond(mt.sample_qc[measure] < cutoff_lower,
                                                             mt.failing_samples_qc.append(f"failing_{measure}"),
                                                             mt.failing_samples_qc))

            # Collect thresholds for each batch
            batch_thresholds[batch][measure] = {'min_thresh': cutoff_lower, 'max_thresh': cutoff_upper}

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
        p = hl.plot.scatter(mt_cols.plot_batch, mt_cols.sample_qc[measure], label=mt_cols[args.batch_col_name],
                            title=f"{measure} values split by batch.")
        save(p)

    ##########################
    # Report failing samples #
    ##########################
    for measure in ['r_ti_tv', 'r_het_hom_var', 'r_insertion_deletion', 'n_singleton']:
        failing_count = mt.aggregate_cols(hl.agg.count_where(mt.failing_samples_qc.contains([f"failing_{measure}"])))

        logging.info(f"Number of samples failing on {measure}: {failing_count}")

    failing_any = mt.aggregate_cols(hl.agg.count_where(hl.len(mt.failing_samples_qc != 0)))
    logging.info(f"Number of samples failing samples QC on any measure: {failing_any}")

    if args.pheno_col is not None:
        cases_failing = mt.aggregate_cols(hl.agg.filter(mt[args.pheno_col] == True,
                                                           hl.agg.count_where(hl.len(mt.failing_samples_qc != 0))))
        controls_failing = mt.aggregate_cols(hl.agg.filter(mt[args.pheno_col] == False,
                                                           hl.agg.count_where(hl.len(mt.failing_samples_qc != 0))))
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
                                                      'batches': batch_set,
                                                      'batch_cohort_name': args.batch_col_name})

    mt_to_annotate = mt_to_annotate.annotate_globals(samples_qc_batch_thresholds=batch_thresholds)

    return mt_to_annotate