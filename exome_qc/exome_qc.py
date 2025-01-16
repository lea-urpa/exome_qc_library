"""
Script for doing exome sequencing data quality control from a Hail matrix table input that has been VEP annotated.

Author: Lea Urpa, August 2020
This pipeline is dedicated to Thao and the Get Down Stay Down, whose album A Man Alive was the mojo for writing
"""
import logging
import os
import hail as hl
from parse_arguments import parse_arguments, check_inputs
import utils
import samples_annotation as sa
import variant_qc as vq
import samples_qc as sq
import variant_annotation as va


def initialize_logger(args):
    """Initializes and configures the logger."""
    datestr = time.strftime("%Y.%m.%d")
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")
    args.log_file = f'exome_qc_{timestr}.txt'

    root = logging.getLogger()
    root.setLevel(logging.INFO if not args.log_debug else logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    file_handler = logging.FileHandler(args.log_file)
    file_handler.setFormatter(formatter)
    root.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    root.addHandler(console_handler)


def verify_inputs(args):
    """Verifies that all input arguments are correct, and that the input file GC bucket matches the dataproc cluster"""
    for folder in [args.log_dir, args.out_dir, args.mt]:
        utils.check_regions(args.region, folder)

    check_inputs(args)


def load_and_annotate_samples(args):
    """Loads and annotates samples based on input arguments."""
    checkpoint_path = os.path.join(args.checkpoint_dir, f"samples_annotated.mt")

    if not utils.check_exists(checkpoint_path) or args.force:
        logging.info(f"Loading matrix table: {args.mt}")
        mt = hl.read_matrix_table(args.mt)
        mt.annotate_globals(original_mt_input={'file': args.mt})
        utils.check_vep(mt)
        utils.check_entry(mt)

        if args.test:
            utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)
            mt = utils.create_test_dataset(mt, args.reference_genome, args.mt, args.out_dir)
            utils.remove_secondary(args.cluster_name, args.region)

        if args.samples_annotation_files:
            logging.info('Annotating samples.')
            for annotation_file in args.samples_annotation_files.split(","):
                mt = sa.annotate_cols_from_file(
                    mt,
                    annotation_file,
                    args.samples_delim,
                    args.samples_col,
                    args.samples_miss
                )

        if args.bam_metadata:
            mt = sa.annotate_cols_from_file(
                mt,
                args.bam_metadata,
                args.bam_delim,
                args.bam_sample_col,
                args.bam_miss
            )

        mt = mt.checkpoint(checkpoint_path, overwrite=True)
    else:
        logging.info("Detected sample-annotated mt exists, skipping samples annotation.")
        mt = hl.read_matrix_table(checkpoint_path)

    return mt


def verify_annotations(args):
    """Confirms necessary annotations exist in the data after annotation."""
    for column_name in args.sample_cols_check:
        col = getattr(args, column_name)
        try:
            test = hl.is_defined(mt[col])
        except Exception as e:
            logging.error(f"Error! Given column annotation {col} does not actually exist after inputting sample "
                          f"annotations.")
            logging.error(e)
            exit(1)

def run_variant_qc(mt, args):
    """Performs variant quality control."""
    qc_path = os.path.join(args.out_dir, "variant_qc.mt")

    if not utils.check_exists(qc_path) or args.force:
        mt = vq.variant_quality_control(
            mt,
            qc_path,
            annotation_prefix="low_pass",
            min_dp=args.min_dp,
            min_gq=args.min_gq,
            max_het_ref_reads=args.max_het_ref_reads,
            min_het_ref_reads=args.min_het_ref_reads,
            call_rate=args.low_pass_min_call_rate,
            p_hwe=args.low_pass_p_hwe,
            snp_qd=args.snp_qd,
            indel_qd=args.indel_qd,
            ab_allowed_dev_het=args.ab_allowed_dev_het,
            count_failing=args.count_failing
        )
    else:
        mt = hl.read_matrix_table(qc_path)

    return mt





def main():
    hl.init()

    args = parse_arguments(sys.argv[1:])
    initialize_logger(args)
    verify_inputs(args)

    mt = load_and_annotate_samples(args)
    mt = run_variant_qc(mt, args)
    mt = run_samples_qc(mt, args)

    mt.write(os.path.join(args.out_dir, "final_qc.mt"), overwrite=True)

if __name__ == "__main__":
    main()


import sys
import time
from bokeh.io import output_file, save

if __name__ == "__main__":





    low_pass_qcd = os.path.join(args.out_dir, f"{stepcount}_{args.out_name}_low_pass_qcd{args.test_str}.mt/")

    #######################
    # low-pass variant QC #
    #######################
    if (not utils.check_exists(low_pass_qcd)) or args.force:
        logging.info("Running low-pass variant QC and genotype QC before samples QC.")

        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)
        mt = hl.read_matrix_table(samples_cleaned)

        mt = vq.variant_quality_control(
            mt, low_pass_qcd, annotation_prefix="low_pass", min_dp=args.min_dp, min_gq=args.min_gq,
            max_het_ref_reads=args.max_het_ref_reads, min_het_ref_reads=args.min_het_ref_reads,
            min_hom_ref_ref_reads=args.min_hom_ref_ref_reads, max_hom_alt_ref_reads=args.max_hom_alt_ref_reads,
            call_rate=args.low_pass_min_call_rate, p_hwe=args.low_pass_p_hwe, snp_qd=args.snp_qd, indel_qd=args.indel_qd,
            ab_allowed_dev_het=args.ab_allowed_dev_het,count_failing=args.count_failing, sex_aware_call_rate=False,
            pheno_col=args.pheno_col, samples_qc=False, force=args.force
        )

        mt = mt.checkpoint(low_pass_qcd, overwrite=True)

        output_file(f"low_pass_mean_het_ab_hist.html")
        ab_hist = mt.aggregate_rows(hl.agg.hist(mt.low_pass_het_ab_stats.mean, 0, 1, 50))
        p = hl.plot.histogram(ab_hist, legend='het ref read ratio', title='Mean het read ratio per var (passing GTs)')
        save(p)

        output_file(f"low_pass_initial_call_rate.html")
        cr_hist_1 = mt.aggregate_rows(hl.agg.hist(mt.low_pass_initial_call_rate, 0, 1, 50))
        p1 = hl.plot.histogram(cr_hist_1, legend='call rate', title="variant call rate, before GT filters")
        save(p1)

        logging.info(f"Writing checkpoint {stepcount}: low pass variant QC")

        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
    else:
        logging.info("Detected low-pass variant QC mt exists, skipping low-pass variant QC.")

    stepcount += 1
    relatedness_calculated = os.path.join(
        args.out_dir, f"{stepcount}_{args.out_name}_relatedness_calculated{args.test_str}.mt/")

    #########################
    # Calculate relatedness #
    #########################
    ld_pruned = os.path.join(args.out_dir, f"{stepcount}-1_{args.out_name}_ld_pruned{args.test_str}.mt/")
    ld_pruned_maffilt = os.path.join(args.out_dir, f"{stepcount}-2_{args.out_name}_maf_filt{args.test_str}.mt/")
    ld_pruned_annot = os.path.join(args.out_dir, f"{stepcount}-2_{args.out_name}_ld_pruned_related{args.test_str}.mt/")

    if (not utils.check_exists(relatedness_calculated)) or args.force:
        logging.info("Calculating relatedness")
        mt = hl.read_matrix_table(low_pass_qcd)
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        ## LD prune and checkpoint ##
        if (not utils.check_exists(ld_pruned)) or args.force:
            logging.debug(f"MAF distribution for variants, before MAF filtering and before filtering failing variants:"
                          f" {mt.low_pass_variant_qc.AF[1].summarize()}")

            # Filter failing samples, variants, and genotypes
            mt_gt_filt = sq.filter_failing(
                mt, ld_pruned, prefix='low_pass', entries=True, variants=False, samples=False, unfilter_entries=True,
                pheno_qc=False, min_dp=args.min_dp, min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
                min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
                max_hom_alt_ref_reads=args.max_hom_alt_ref_reads, force=args.force
            )

            mt_filtered = mt.filter_rows(hl.len(mt.low_pass_failing_variant_qc) == 0)
            logging.debug(f"Number of variants in low-pass QCd matrix table: {mt_filtered.count_rows()}")
            logging.debug(f"MAF distribution for variants, before MAF filtering but after filtering failing variants:"
                          f" {mt_filtered.low_pass_variant_qc.AF[1].summarize()}")

            # Filter out low MAF variants
            if (not utils.check_exists(ld_pruned_maffilt)) or args.force:
                mt_maffilt = vq.maf_filter(mt_filtered, args.ind_maf, "low_pass_variant_qc")
                mt_maffilt = mt_maffilt.checkpoint(ld_pruned_maffilt, overwrite=True)
                logging.debug(f"Number of variants in MAF filtered {args.ind_maf} matrix table: {mt_maffilt.count_rows()}")
            else:
                mt_maffilt = hl.read_matrix_table(ld_pruned_maffilt)
                logging.debug(f"Number of variants in MAF filtered {args.ind_maf} matrix table: {mt_maffilt.count_rows()}")

            # LD prune if row count >80k
            mt_ldpruned = vq.downsample_variants(
                mt_maffilt, 80000, ld_pruned, r2=args.r2, bp_window_size=args.bp_window_size, ld_prune=True)

            logging.info(f"Writing checkpoint {stepcount}-1: LD pruned dataset")
            logging.debug(f"Count of variants in downsampled matrix table: {mt_ldpruned.count_rows()}")
            mt_ldpruned = mt_ldpruned.checkpoint(ld_pruned, overwrite=True)
        else:
            logging.info("Detected LD pruned dataset written, loading that.")
            logging.debug(f"Count of variants in downsampled matrix table: {mt_ldpruned.count_rows()}")
            mt_ldpruned = hl.read_matrix_table(ld_pruned)

        ## Calculate relatedness with King ##
        if args.reference_genome == "GRCh38":
            autosomes = ["chr" + str(i) for i in range(1, 23)]
        else:
            autosomes = [str(i) for i in range(1, 23)]

        mt_autosomes = mt_ldpruned.filter_rows(hl.literal(autosomes).contains(mt_ldpruned.locus.contig))

        related_to_remove, related_info_ht = sq.king_relatedness(
            mt_autosomes, relatedness_calculated, kinship_threshold=args.kinship_threshold, pheno_col=args.pheno_col,
            force=args.force, cluster_name=args.cluster_name, num_secondary_workers=args.num_secondary_workers,
            region=args.region)

        mt = mt.annotate_cols(
            related_to_remove=hl.if_else(hl.literal(related_to_remove).contains(mt.s), True, False),
            related_graph_id=related_info_ht[mt.s].related_graph_id,
            related_num_connections=related_info_ht[mt.s].related_num_connections
        )

        mt = mt.annotate_cols(related_num_connections=hl.or_else(mt.related_num_connections, 0))

        mt_ldpruned = mt_ldpruned.annotate_cols(
            related_to_remove=hl.if_else(hl.literal(related_to_remove).contains(mt_ldpruned.s), True, False),
            related_graph_id=related_info_ht[mt_ldpruned.s].related_graph_id,
            related_num_connections=related_info_ht[mt_ldpruned.s].related_num_connections)

        mt_ldpruned = mt_ldpruned.annotate_cols(related_num_connections=hl.or_else(mt_ldpruned.related_num_connections, 0))

        logging.info(f"Writing checkpoint {stepcount}: relatedness annotated")
        mt = mt.checkpoint(relatedness_calculated, overwrite=True)
        mt_ldpruned = mt_ldpruned.checkpoint(ld_pruned_annot, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

    else:
        logging.info("Detected mt with relatives annotated exists, skipping relatedness calculation.")

    stepcount += 1
    pop_outliers_found = os.path.join(
        args.out_dir, f"{stepcount}_{args.out_name}_pop_outliers_found{args.test_str}.mt/")

    ############################
    # Find population outliers #
    ############################
    ld_pruned_maf = os.path.join(args.out_dir, f"{stepcount}-1_{args.out_name}_ld_pruned_maf_0.05{args.test_str}.mt/")
    ld_pruned_popannot = os.path.join(args.out_dir,
                                      f"{stepcount}-2_{args.out_name}_ld_pruned_popoutliers{args.test_str}.mt/")

    if (not utils.check_exists(pop_outliers_found)) or args.force:
        logging.info("Finding population outliers")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        mt = hl.read_matrix_table(relatedness_calculated)
        mt_ldpruned = hl.read_matrix_table(ld_pruned_annot)

        maf_stats = mt_ldpruned.aggregate_rows(hl.agg.stats(mt_ldpruned.low_pass_variant_qc.AF[1]))

        if maf_stats.min < 0.05:
            logging.info("Further excluding variants with MAF < 0.05 to calculate principal components.")
            mt_ldpruned = mt_ldpruned.filter_rows(mt_ldpruned.variant_qc.AF[1] >= 0.05, keep=True)
            mt_ldpruned = mt_ldpruned.checkpoint(ld_pruned_maf, overwrite=True)

        pop_outliers = sq.find_pop_outliers(
            mt_ldpruned, pop_outliers_found, pop_sd_threshold=args.pop_sd_threshold,
            plots=args.pca_plots, max_iter=args.max_iter, reference_genome=args.reference_genome,
            pca_plot_annotations=args.pca_plot_annotations)

        mt = mt.annotate_cols(pop_outlier_sample=hl.if_else(hl.literal(pop_outliers).contains(mt.s), True, False))
        mt_ldpruned = mt_ldpruned.annotate_cols(
            pop_outlier_sample=hl.if_else(hl.literal(pop_outliers).contains(mt_ldpruned.s), True, False))

        logging.info(f"Writing checkpoint {stepcount}: population outliers annotated")
        mt = mt.checkpoint(pop_outliers_found, overwrite=True)
        mt_ldpruned = mt_ldpruned.checkpoint(ld_pruned_popannot, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

    else:
        logging.info("Detected mt with population outliers annotated exists, skipping finding pop outliers.")

    stepcount += 1
    sex_imputed = os.path.join(args.out_dir, f"{stepcount}_{args.out_name}_sex_imputed{args.test_str}.mt/")

    ####################################
    # Annotate variants and impute sex #
    ####################################
    filtered_nohwe = os.path.join(args.out_dir, f"{stepcount}-1_{args.out_name}_filtered_except_hwe{args.test_str}.mt/")
    filtered_annot = os.path.join(args.out_dir, f"{stepcount}-2_{args.out_name}_filtered_annotated{args.test_str}.mt/")
    unfilt_annot = os.path.join(args.out_dir, f"{stepcount}-3_{args.out_name}_annotated_tmp")

    if (not utils.check_exists(sex_imputed)) or args.force:
        logging.info("Annotating variants and imputing sex")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        mt = hl.read_matrix_table(pop_outliers_found)

        # Annotate variants
        mt = va.annotate_variants(mt)

        if (not utils.check_exists(filtered_nohwe)) or args.force:
            # Filter out failing variants, genotypes, rare variants
            mt_gt_filt= sq.filter_failing(
                mt, sex_imputed, prefix='low_pass', variants=False, entries=True, samples=False,
                unfilter_entries=False, pheno_qc=False, min_dp=args.min_dp,
                min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
                min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
                max_hom_alt_ref_reads=args.max_hom_alt_ref_reads, force=args.force
            )

            mt_filtered = mt_gt_filt.filter_rows(
                (mt_gt_filt.low_pass_failing_variant_qc == ["failing_hwe"]) |
                (hl.len(mt_gt_filt.low_pass_failing_variant_qc) == 0), keep=True
            )

            mt_filtered = mt_filtered.checkpoint(filtered_nohwe, overwrite=True)
        else:
            mt_filtered = hl.read_matrix_table(filtered_nohwe)

        # Impute sex
        imputed_sex = sq.impute_sex_plot(mt_filtered, female_threshold=args.female_threshold,
                                         male_threshold=args.male_threshold, aaf_threshold=0.05)

        # Annotate unfiltered + variant/GT filtered mt with imputed sex values
        mt_filtered = mt_filtered.annotate_cols(is_female_imputed=imputed_sex[mt_filtered.s].is_female)
        mt = mt.annotate_cols(is_female_imputed=imputed_sex[mt.s].is_female, f_stat=imputed_sex[mt.s].f_stat)
        mt = mt.annotate_globals(
            sex_imputation_thresholds={'female_threshold': args.female_threshold,'male_threshold': args.male_threshold})

        # Calculate sex-aware sample annotations with variant filtered mt
        mt_filtered = sa.sex_aware_sample_annotations(mt_filtered)
        mt_filtered = mt_filtered.checkpoint(filtered_annot, overwrite=True)

        mt = mt.annotate_cols(sexaware_sample_call_rate=mt_filtered.cols()[mt.s].sexaware_sample_call_rate)

        logging.info(f"Writing checkpoint {stepcount}: variants annotated and sex imputed")
        mt = mt.checkpoint(sex_imputed, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
    else:
        logging.info("Detected variant annotation and sex imputation completed, skipping this step.")

    stepcount += 1
    samples_qcd = os.path.join(args.out_dir, f"{stepcount}_{args.out_name}_samples_qcd{args.test_str}.mt/")

    ##############
    # Samples QC #
    ##############
    if (not utils.check_exists(samples_qcd)) or args.force:
        logging.info("Running sample QC")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        mt = hl.read_matrix_table(sex_imputed)

        # Filter failing variants and genotypes, and pop outlier samples
        logging.info("Filtering out failing genotypes and variants.")
        mt_filtered = sq.filter_failing(
            mt, samples_qcd, prefix='low_pass', entries=True, variants=True, samples=False, unfilter_entries=False,
            pheno_qc=False, min_dp=args.min_dp, min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
            min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
            max_hom_alt_ref_reads=args.max_hom_alt_ref_reads, force=args.force
        )
        utils.remove_secondary(args.cluster_name, args.region)

        # Run samples QC
        mt = sq.samples_qc(
            mt_filtered, mt, samples_qcd, count_failing=args.count_failing, sample_call_rate=args.sample_call_rate,
            chimeras_col=args.chimeras_col, chimeras_max=args.chimeras_max, contamination_col=args.contamination_col,
            contamination_max=args.contamination_max, batch_col_name=args.batch_col_name,
            sampleqc_sd_threshold=args.sampleqc_sd_threshold, pheno_col=args.pheno_col
        )

        logging.info(f"Writing checkpoint {stepcount}: sample QC")
        mt = mt.checkpoint(samples_qcd, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
    else:
        logging.info("Detected samples QC completed, skipping this step.")

    stepcount += 1
    variant_qcd = os.path.join(args.out_dir, f"{stepcount}_{args.out_name}_variant_qcd{args.test_str}.mt/")

    ##############
    # Variant QC #
    ##############
    if (not utils.check_exists(variant_qcd)) or args.force:
        logging.info("Running final variant QC")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        mt = hl.read_matrix_table(samples_qcd)

        # Run variant QC
        mt = vq.variant_quality_control(
            mt, variant_qcd, annotation_prefix="final", min_dp=args.min_dp, min_gq=args.min_gq,
            max_het_ref_reads=args.max_het_ref_reads, min_het_ref_reads=args.min_het_ref_reads,
            min_hom_ref_ref_reads=args.min_hom_ref_ref_reads, max_hom_alt_ref_reads=args.max_hom_alt_ref_reads,
            call_rate=args.final_min_call_rate, p_hwe=args.final_p_hwe, snp_qd=args.snp_qd, indel_qd=args.indel_qd,
            ab_allowed_dev_het=args.ab_allowed_dev_het,
            count_failing=args.count_failing, sex_aware_call_rate=True, pheno_col=args.pheno_col,
            samples_qc=True, force=args.force
        )

        # Run case status-specific variant QC
        if args.pheno_col is not None:
            logging.info("Running case-status specific variant QC")
            mt = vq.find_variants_failing_by_pheno(mt, ab_allowed_dev_het=args.ab_allowed_dev_het,
                                                   pheno_call_rate=args.pheno_call_rate, prefix = "final")
        else:
            logging.info("Phenotype column not given, skipping filtering variants by phenotype.")

        logging.info(f"Writing checkpoint {stepcount}: variant QC")
        mt = mt.checkpoint(variant_qcd, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
    else:
        logging.info("Detected that final variant QC has been run. Skipping this step.")

    stepcount += 1
    pcs_calculated = os.path.join(args.out_dir, f"{stepcount}_{args.out_name}_final_with_PCs{args.test_str}.mt/")

    #######################
    # Calculate final PCs #
    #######################
    final_filtered = os.path.join(args.out_dir, f"{stepcount}-1_{args.out_name}_filtered{args.test_str}.mt/")
    final_maffilt = os.path.join(args.out_dir, f"{stepcount}-2_{args.out_name}_maf_filt{args.test_str}.mt/")
    final_ldpruned = os.path.join(args.out_dir, f"{stepcount}-3_{args.out_name}_ldpruned{args.test_str}.mt/")

    if (not utils.check_exists(pcs_calculated)) or args.force:
        logging.info("Calculating final PCs")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        mt = hl.read_matrix_table(variant_qcd)

        # Filter out failing samples, variants, genotypes for PC calculations
        if (not utils.check_exists(final_filtered)) or args.force:
            mt_filtered = sq.filter_failing(
                mt, pcs_calculated, prefix='final', entries=True, variants=True, samples=True, unfilter_entries=True,
                pheno_qc=False, min_dp=args.min_dp, min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
                min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
                max_hom_alt_ref_reads=args.max_hom_alt_ref_reads, force=args.force, pop_outliers=True
            )
            mt_filtered = mt_filtered.checkpoint(final_filtered, overwrite=True)

        else:
            logging.info("Detected final failing sample, variant, and genotype mt exists. Loading that.")
            mt_filtered = hl.read_matrix_table(final_filtered)

        # MAF filter
        if (not utils.check_exists(final_maffilt)) or args.force:
            logging.info("Filtering to common variants and LD pruning dataset.")
            mt_maffilt = vq.maf_filter(mt_filtered, args.ind_maf, "final_variant_qc")
            mt_maffilt = mt_maffilt.checkpoint(final_maffilt, overwrite=True)
        else:
            logging.info("Detected final MAF filtered mt exists. Loading that.")
            mt_maffilt = hl.read_matrix_table(final_maffilt)

        # LD prune
        if (not utils.check_exists(final_ldpruned)) or args.force:
            logging.info('LD pruning final dataset for PC calculation')
            mt_ldpruned = vq.downsample_variants(
                mt_maffilt, 80000, final_ldpruned, r2=args.r2, bp_window_size=args.bp_window_size, ld_prune=True)
            mt_ldpruned = mt_ldpruned.checkpoint(final_ldpruned, overwrite=True)
        else:
            logging.info("Detected final LDpruned mt exists. Loading that.")
            mt_ldpruned = hl.read_matrix_table(final_ldpruned)

        # Calculate PCs and project to relatives, annotate to main mt, plot
        scores, related_scores = sq.project_pcs_relateds(mt_ldpruned, pcs_calculated, args.pc_num, args.reference_genome)

        mt = mt.annotate_cols(**{'pc' + str(k + 1): scores[mt.s].scores[k]
                                 for k in range(args.pc_num)})
        mt = mt.annotate_cols(**{'pc' + str(k + 1): hl.or_else(mt['pc' + str(k + 1)], related_scores[mt.s].scores[k])
                                 for k in range(args.pc_num)})

        if args.pca_plot_annotations is not None:
            try:
                pca_annotations = args.pca_plot_annotations.strip().split(",")
                label_dict = {i: mt[i] for i in pca_annotations}

                output_file(f"{datestr}_final_pcs_plot.html")
                p = hl.plot.scatter(mt.pc1, mt.pc2, label=label_dict, title="Final PCs", collect_all=True)
                save(p)
            except Exception as e:
                logging.error(f"Error! Creating PCA plots with labels failed. Are the label categories you provided"
                              f" really in the data? labels provided: {args.pca_plot_annotations}. Plotting without "
                              f"labels")
                logging.error(e)
                output_file(f"{datestr}_final_pcs_plot.html")
                p = hl.plot.scatter(mt.pc1, mt.pc2, title="Final principal components", collect_all=True)
                save(p)
        else:
            output_file(f"{datestr}_final_pcs_plot.html")
            p = hl.plot.scatter(mt.pc1, mt.pc2, title="Final principal components", collect_all=True)
            save(p)

        mt = mt.checkpoint(pcs_calculated, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

    else:
        logging.info("Detected that final PCs have been calculated. Skipping this step.")

    stepcount += 1
    variant_annotated = os.path.join(args.out_dir, f"{stepcount}_{args.out_name}_variants_annotated{args.test_str}.mt/")

    ##############################
    # Annotate with CADD, Gnomad #
    ##############################
    if (not utils.check_exists(variant_annotated)) or args.force:
        mt = hl.read_matrix_table(pcs_calculated)

        if args.mpc_ht is not None:
            logging.info("Annotating variants with MPC info.")
            mt = va.annotate_variants_mpc(mt, args.mpc_ht)
        if args.cadd_ht is not None:
            logging.info("Annotating variants with CADD info.")
            mt = va.annotate_variants_cadd(mt, args.cadd_ht)
        if args.gnomad_ht is not None:
            logging.info("Annotating variants with Gnomad.")
            mt = va.annotate_variants_gnomad(mt, args.gnomad_ht)
        if args.gnomad_mismatch_ht is not None:
            logging.info("Annotating variants with gnomad mismatch info.")
            mt = va.annotate_variants_gnomad_mismatch(mt, args.gnomad_mismatch_ht)

        mt = mt.checkpoint(variant_annotated, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
    else:
        mt = hl.read_matrix_table(variant_annotated)
        logging.info("Detected final annotated mt exists. Skipping this step.")

    # Export rows and columns
    mtcols = mt.cols()
    mtcols = mtcols.flatten()
    mtcols.export(os.path.join(args.out_dir, args.out_name + '_final_dataset_cols.tsv.gz'))

    mtrows = mt.rows()
    mtrows = mtrows.flatten()
    mtrows = mtrows.key_by().drop("vep.input")
    mtrows.export(os.path.join(args.out_dir, args.out_name + '_final_dataset_rows.tsv.gz'))


    # Send logs and finish-up notice
    logging.info('Pipeline ran successfully! Copying logs and shutting down cluster in 10 minutes.')
    utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

