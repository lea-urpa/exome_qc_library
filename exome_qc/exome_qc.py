"""
Script for doing exome sequencing data quality control from a Hail matrix table input that has been VEP annotated.

Author: Lea Urpa, August 2020
"""
import sys
import os
import logging
import time
import hail as hl
from bokeh.io import output_file, save
from parse_arguments import parse_arguments, check_inputs
import utils
import samples_annotation as sa
import variant_qc as vq
import samples_qc as sq
import variant_annotation as va

if __name__ == "__main__":
    ###################################################################
    # Initialize Hail and import scripts, configure logger and inputs #
    ###################################################################
    hl.init()

    args = parse_arguments(sys.argv[1:])
    check_inputs(args)

    ## Configure logger ##
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day
    args.log_file = 'exome-qc_' + timestr + '.txt'

    root = logging.getLogger() # creates logger
    root.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Add file handler
    fh = logging.FileHandler(args.log_file)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add streaming handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    ## Configure inputs ##
    args.checkpoint_folder = os.path.join(args.out_dir, "checkpoint_mts/")
    args.plot_folder = os.path.join(args.out_dir, "plots")
    stepcount = 1
    if args.test:
        args.test_str = "_test"
    else:
        args.test_str = ""

    ##################################
    # Load data and annotate samples #
    ##################################
    samples_annotated = os.path.join(args.out_dir, f"{stepcount}_{args.out_name}_samples_annotated{args.test_str}.mt/")

    if (not utils.check_exists(samples_annotated)) or args.force:

        ## Load data ##
        logging.info(f"Loading matrix table: {args.mt}")
        mt = hl.read_matrix_table(args.mt)
        mt.annotate_globals(original_mt_input={'file': args.mt, 'date': datestr})

        utils.check_vep(mt)


        if args.test:
            utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)
            mt = utils.create_test_dataset(mt, args.reference_genome, args.mt, args.out_dir)
            utils.remove_secondary(args.cluster_name, args.region)

        ## Annotate samples ##
        logging.info('Annotating samples.')

        # Annotate with optional samples annotation files
        if args.samples_annotation_files is not None:
            annotation_files = args.samples_annotation_files.strip().split(",")

            for file in annotation_files:
                mt = sa.annotate_cols_from_file(mt, file, args.samples_delim, args.samples_col, args.samples_miss)

        # Annotate with bam metadata
        if args.bam_metadata is not None:
            mt = sa.annotate_cols_from_file(mt, args.bam_metadata, args.bam_delim, args.bam_sample_col, args.bam_miss)

        # Check columns exist
        for colname in args.sample_cols_check:
            col = getattr(args, colname)
            try:
                test = hl.is_defined(mt[col])
            except Exception as e:
                logging.error(f"Error! Given column annotation {col} does not actually exist after inputting sample "
                              f"annotations.")
                logging.error(e)
                exit(1)

        logging.info(f"Writing checkpoint {stepcount}: annotating samples")
        mt = mt.checkpoint(samples_annotated, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

    else:
        logging.info("Detected sample-annotated mt exists, skipping samples annotation.")

    samples_removed = os.path.join(args.out_dir, f"{stepcount}-1_{args.out_name}_samples_removed{args.test_str}.mt/")

    ##################
    # Remove samples #
    ##################
    if (not utils.check_exists(samples_removed) or args.force) and (
            (args.sample_removal_strings is not None) or (args.sample_removal_list is not None)):

        logging.info("Removing indicated samples.")
        mt = hl.read_matrix_table(samples_annotated)

        samples_start = mt.count_cols()
        logging.info(f"Initial sample count: {samples_start}")

        # Filter out samples from arbitrary list uploaded with args
        if args.sample_removal_list is not None:
            rm_list = hl.import_table(args.sample_removal_list, no_header=True)
            rm_list = rm_list.annotate(s=rm_list.f0)
            rm_list = rm_list.key_by("s")
            mt = mt.anti_join_cols(rm_list)
            list_filtered = mt.count_cols()
            logging.info(f"Sample count after filtering by input list: {list_filtered}")

        # Filter out samples that have sample names containing a particular string
        if args.sample_removal_strings is not None:
            removal_strings = args.sample_removal_strings.strip().split(",")
            for key_string in removal_strings:
                mt = mt.filter_cols(mt.s.contains(key_string), keep=False)
            string_filtered = mt.count_cols()
            logging.info(f"Sample count after filtering by strings: {string_filtered}")

        logging.info(f"Writing checkpoint {stepcount}-1: sample removal")
        mt = mt.checkpoint( samples_removed, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

        samples_cleaned = samples_removed
    else:
        samples_cleaned = samples_annotated

    stepcount += 1
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

        logging.info(f"Writing checkpoint {stepcount}: low pass variant QC")
        mt = mt.checkpoint(low_pass_qcd, overwrite=True)
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
    ld_pruned_annot = os.path.join(args.out_dir, f"{stepcount}-2_{args.out_name}_ld_pruned_related{args.test_str}.mt/")

    if (not utils.check_exists(relatedness_calculated)) or args.force:
        logging.info("Calculating relatedness")
        mt = hl.read_matrix_table(low_pass_qcd)

        ## LD prune and checkpoint ##
        if not utils.check_exists(ld_pruned):
            utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

            # Filter failing samples, variants, and genotypes
            mt_filtered = sq.filter_failing(
                mt, ld_pruned, prefix='low_pass', entries=True, variants=True, samples=False, unfilter_entries=True,
                pheno_qc=False, min_dp=args.min_dp, min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
                min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
                max_hom_alt_ref_reads=args.max_hom_alt_ref_reads
            )

            # Filter out low MAF variants
            mt_maffilt = vq.maf_filter(mt_filtered, args.ind_maf, filter_ac0_after_pruning=True)

            # LD prune if row count >80k
            utils.remove_secondary(args.cluster_name, args.region)
            mt_ldpruned = vq.downsample_variants(mt_maffilt, 80000, r2=args.r2, bp_window_size=args.bp_window_size,
                                                 ld_prune=True)

            logging.info(f"Writing checkpoint {stepcount}-1: LD pruned dataset")
            mt_ldpruned = mt_ldpruned.checkpoint(ld_pruned, overwrite=True)
        else:
            logging.info("Detected LD pruned dataset written, loading that.")
            mt_ldpruned = hl.read_matrix_table(ld_pruned)

        ## Calculate relatedness with King ##
        if not utils.check_exists(relatedness_calculated):

            if args.reference_genome is "GRCh38":
                autosomes = ["chr" + str(i) for i in range(1, 23)]
            else:
                autosomes = [str(i) for i in range(1, 23)]

            mt_autosomes = mt_ldpruned.filter_rows(hl.literal(autosomes).contains(mt_ldpruned.locus.contig))

            related_to_remove = sq.king_relatedness(
                mt_autosomes, relatedness_calculated, kinship_threshold=args.kinship_threshold, pheno_col=args.pheno_col,
                force=args.force, cluster_name=args.cluster_name, num_secondary_workers=args.num_secondary_workers,
                region=args.region, reference_genome=args.reference_genome)

            mt = mt.annotate_cols(
                related_to_remove=hl.if_else(hl.literal(related_to_remove).contains(mt.s), True, False))

            mt_ldpruned = mt_ldpruned.annotate_cols(
                related_to_remove=hl.if_else(hl.literal(related_to_remove).contains(mt_ldpruned.s), True, False))

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
    ld_pruned_popannot = os.path.join(args.out_dir,
                                      f"{stepcount}-1_{args.out_name}_ld_pruned_popoutliers{args.test_str}.mt/")

    if (not utils.check_exists(pop_outliers_found)) or args.force:
        logging.info("Finding population outliers")

        mt = hl.read_matrix_table(relatedness_calculated)
        mt_ldpruned = hl.read_matrix_table(ld_pruned_annot)

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
    filtered_annot = os.path.join(args.out_dir, f"{stepcount}-2_{args.out_name}_filtered_sex_annotated{args.test_str}.mt/")

    if (not utils.check_exists(sex_imputed)) or args.force:
        logging.info("Annotating variants and imputing sex")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        mt = hl.read_matrix_table(pop_outliers_found)

        # Annotate variants
        mt = va.annotate_variants(mt)

        if (not utils.check_exists(filtered_nohwe)) or args.force:
            # Filter out failing variants, genotypes, rare variants
            mt_filtered = sq.filter_failing(
                mt, sex_imputed, prefix='low_pass', variants=True, entries=True, samples=False,
                unfilter_entries=False, skip_hwe_fails=True, pheno_qc=False, min_dp=args.min_dp,
                min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
                min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
                max_hom_alt_ref_reads=args.max_hom_alt_ref_reads
            )

            mt_filtered = mt_filtered.checkpoint(filtered_nohwe, overwrite=True)
        else:
            mt_filtered = hl.read_matrix_table(filtered_nohwe)

        # Impute sex
        imputed_sex = sq.impute_sex_plot(mt_filtered, female_threshold=args.female_threshold,
                                         male_threshold=args.male_threshold, aaf_threshold=0.05)

        # Annotate unfiltered + filtered mt with imputed sex values
        mt_filtered = mt_filtered.annotate_cols(is_female_imputed=imputed_sex[mt_filtered.s].is_female)
        mt = mt.annotate_cols(is_female_imputed=imputed_sex[mt.s].is_female, f_stat=imputed_sex[mt.s].f_stat)
        mt = mt.annotate_globals(
            sex_imputation_thresholds={'female_threshold': args.female_threshold,'male_threshold': args.male_threshold})

        # Annotate sex-aware variant annotations (gt filt only to have for all
        gt_filt_fn = sex_imputed.rstrip("/").replace(".mt", "") + "_GT_filtered.mt/"
        mt_gt_filt = hl.read_matrix_table(gt_filt_fn)
        mt_gt_filt = mt_gt_filt.annotate_cols(is_female_imputed=imputed_sex[mt_gt_filt.s].is_female)

        mt_filtered, annotations_to_transfer = va.sex_aware_variant_annotations(mt_gt_filt, pheno_col=args.pheno_col)
        mt_filtered = sa.sex_aware_sample_annotations(mt_filtered)

        mt_filtered = mt_filtered.checkpoint(filtered_annot, overwrite=True)

        for annotation in annotations_to_transfer:
            mt = mt.annotate_rows(**{annotation: mt_filtered.rows()[mt.row_key][annotation]})
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
        utils.remove_secondary(args.cluster_name, args.region)

        mt = hl.read_matrix_table(sex_imputed)

        # Filter failing variants and genotypes, and pop outlier samples
        logging.info("Filtering out failing genotypes and variants.")
        mt_filtered = sq.filter_failing(
            mt, ld_pruned, prefix='low_pass', entries=True, variants=True, samples=False, unfilter_entries=False,
            pheno_qc=False, min_dp=args.min_dp, min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
            min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
            max_hom_alt_ref_reads=args.max_hom_alt_ref_reads
        )

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
        logging.info("Running variant QC")
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
                                                   pheno_call_rate=args.pheno_call_rate)
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
                max_hom_alt_ref_reads=args.max_hom_alt_ref_reads
            )
            mt_filtered = mt_filtered.checkpoint(final_filtered, overwrite=True)

        else:
            logging.info("Detected final failing sample, variant, and genotype mt exists. Loading that.")
            mt_filtered = hl.read_matrix_table(final_filtered)

        # MAF filter
        if (not utils.check_exists(final_maffilt)) or args.force:
            logging.info("Filtering to common variants and LD pruning dataset.")
            mt_maffilt = vq.maf_filter(mt_filtered, 0.05)
            mt_maffilt = mt_maffilt.checkpoint(final_maffilt)
        else:
            logging.info("Detected final MAF filtered mt exists. Loading that.")
            mt_maffilt = hl.read_matrix_table(final_maffilt)

        # LD prune
        if (not utils.check_exists(final_ldpruned)) or args.force:
            logging.info('LD pruning final dataset for PC calculation')
            utils.remove_secondary(args.cluster_name, args.region)
            mt_ldpruned = vq.downsample_variants(mt_maffilt, 80000, r2=args.r2, bp_window_size=args.bp_window_size,
                                                 ld_prune=True)
            mt_ldpruned = mt_ldpruned.checkpoint(final_ldpruned)
        else:
            logging.info("Detected final LDpruned mt exists. Loading that.")
            mt_ldpruned = hl.read_matrix_table(final_ldpruned)

        # Calculate PCs and project to relatives, plot
        mt = sq.project_pcs_relateds(mt_ldpruned, args.pc_num)

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

        # Export rows and columns
        mtcols = mt.cols()
        mtcols = mtcols.flatten()
        mtcols.export(os.path.join(args.out_dir, args.out_name + '_final_dataset_cols.tsv'))

        mtrows = mt.rows()
        mtrows = mtrows.flatten()
        mtrows = mtrows.key_by().drop("vep.input")
        mtrows.export(os.path.join(args.out_dir, args.out_name + '_final_dataset_rows.tsv'))

    else:
        logging.info("Detected final annotated mt exists. Skipping this step.")


    # Send logs and finish-up notice
    logging.info('Pipeline ran successfully! Copying logs and shutting down cluster in 10 minutes.')
    utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
