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

            logging.info('Test flag given, filtering to chrom 22 and chrom X.')
            if args.reference_genome == "GRCh38":
                chrom_codes = hl.array(["chr22", "chrX"])
            else:
                chrom_codes = hl.array(["22", "X"])

            mt = mt.filter_rows(chrom_codes.contains(mt.locus.contig))

            if args.mt.endswith(".mt/"):
                test_mt = (args.mt.replace(".mt/", "") + "_test.mt").split("/")[-1]
            else:
                test_mt = (args.mt[:-1] + "_test.mt").split("/")[-1]
            mt = mt.checkpoint(os.path.join(args.out_dir, test_mt), overwrite=True)

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
            pheno_col=args.pheno_col, samples_qc=False
        )

        logging.info(f"Writing checkpoint {stepcount}: low pass variant QC")
        mt = mt.checkpoint(low_pass_qcd, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
    else:
        logging.info("Detected low-pass variant QC mt exists, skipping low-pass variant QC.")

    stepcount += 1
    relatedness_calculated = os.path.join(
        args.out_dir, f"{stepcount}_{args.out_name}_relatedness_calculated{args.test_str}.mt/")

    logging.info("Tested up to low-pass variant QC, exiting now.")
    exit(0)

    #########################
    # Calculate relatedness #
    #########################
    ld_pruned = os.path.join(args.out_dir, f"{stepcount}-1_{args.out_name}_ld_pruned{args.test_str}.mt/")
    logging.info("This part of the pipeline not implemented yet, exiting now!")
    exit(0)

    if (not utils.check_exists(relatedness_calculated)) or args.force:
        logging.info("Calculating relatedness")
        mt = hl.read_matrix_table(low_pass_qcd)

        ## LD prune and checkpoint ##
        if not utils.check_exists(ld_pruned):
            utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

            # Filter failing samples, variants, and genotypes
            mt_filtered = sq.filter_failing(mt, args, mode='low_pass', unfilter_entries=True, checkpoint=False)

            # Filter out low MAF variants
            mt_maffilt = vq.maf_filter(mt_filtered, args.ind_maf, filter_ac0_after_pruning=True)

            # LD prune if row count >80k
            utils.remove_secondary(args.cluster_name, args.region)
            mt_ldpruned = vq.downsample_variants(mt_maffilt, 80000)

            logging.info(f"Writing checkpoint {stepcount}-1: LD pruned dataset")
            mt_ldpruned = mt_ldpruned.checkpoint(ld_pruned, overwrite=True)
        else:
            logging.info("Detected LD pruned dataset written, loading that.")
            mt_ldpruned = hl.read_matrix_table(ld_pruned)

        ## Calculate relatedness with King ##
        if not utils.check_exists(relatedness_calculated):
            # TODO find out if I can add secondary workers here- does it trigger a shuffle?
            relatedness = hl.king_or_whatever(mt_ldpruned) # TODO figure out how to implement this

            # Annotate with relatedness
            # TODO implement this

            logging.info(f"Writing checkpoint {stepcount}: relatedness annotated")
            mt = mt.checkpoint(relatedness_calculated, overwrite=True)
            utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

    stepcount += 1
    pop_outliers_found = os.path.join(
        args.out_dir, f"{stepcount}_{args.out_name}_relatedness_calculated{args.test_str}.mt/")

    ############################
    # Find population outliers #
    ############################
    if (not utils.check_exists(pop_outliers_found)) or args.force:
        logging.info("Finding population outliers")
        # TODO remove secondary workers if they possibly have been added for King relatedness step

        mt = hl.read_matrix_table(relatedness_calculated)
        mt_ldpruned = hl.read_matrix_table(ld_pruned)

        mt = sq.find_pop_outliers(mt_ldpruned, mt_to_annotate=mt, args=args)

        logging.info(f"Writing checkpoint {stepcount}: population outliers annotated")
        mt = mt.checkpoint(pop_outliers_found, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

    stepcount += 1
    sex_imputed = os.path.join(args.out_dir, f"{stepcount}_{args.out_name}_sex_imputed{args.test_str}.mt/")

    ####################################
    # Annotate variants and impute sex #
    ####################################
    if (not utils.check_exists(sex_imputed)) or args.force:
        logging.info("Annotating variants and imputing sex")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        mt = hl.read_matrix_table(pop_outliers_found)

        # Annotate variants
        mt = va.annotate_variants(mt)

        # Filter out failing samples, rare variants
        mt_filtered = sq.filter_failing(mt, args, mode='low_pass', unfilter_entries=False)
        # TODO fix this so it doesn't take out hardy-weinberg failing variants!
        mt_maf = vq.maf_filter(mt_filtered, 0.05)

        # Impute sex
        #TODO check if this triggers a shuffle
        mt_maf, imputed_sex, mt = sq.impute_sex_plot(mt_maf, mt_to_annotate=mt, args=args)

        # Annotate sex-aware variant annotations
        mt = va.sex_aware_variant_annotations(mt_filtered, mt_to_annotate=mt, args=args)
        mt = sa.sex_aware_sample_annotations(mt_filtered, mt_to_annotate=mt, args=args)
        # TODO check if this triggers a shuffle

        logging.info(f"Writing checkpoint {stepcount}: variants annotated and sex imputed")
        mt = mt.checkpoint(sex_imputed, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

    stepcount += 1
    samples_qcd = os.path.join(args.out_dir, f"{stepcount}_{args.out_name}_samples_qcd{args.test_str}.mt/")

    ##############
    # Samples QC #
    ##############
    if (not utils.check_exists(samples_qcd)) or args.force:
        logging.info("Running sample QC")
        utils.remove_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        mt = hl.read_matrix_table(sex_imputed)

        # Filter failing variants and genotypes, and pop outlier samples
        mt_filtered = sq.filter_failing(mt, args, mode='low_pass', unfilter_entries=False)

        # Run samples QC
        mt = sq.samples_qc(mt_filtered, mt_to_annotate=mt, args=args)

        logging.info(f"Writing checkpoint {stepcount}: sample QC")
        mt = mt.checkpoint(samples_qcd, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

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
        mt = vq.find_failing_variants(mt, args, mode='final')
        # TODO give following parameters
        # Test that failing samples QC and pop outliers detected? Maybe?
        # p_hwe = args.final_p_hwe
        #call_rate = args.final_min_call_rate
        #annotation_name = "variant_qc_thresholds_final"
        #sex_aware_call_rate = "True"
        #varqc_name = 'final_no_failing_samples_varqc'

        # Run case status-specific variant QC
        if args.pheno_col is not None:
            logging.info("Running case-status specific variant QC")
            mt = vq.find_variants_failing_by_pheno(mt, args)
        else:
            logging.info("Phenotype column not given, skipping filtering variants by phenotype.")

        logging.info(f"Writing checkpoint {stepcount}: variant QC")
        mt = mt.checkpoint(variant_qcd, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

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
            mt_filtered = sq.filter_failing(mt, args, mode="final", unfilter_entries=True, pheno_qc=True)
            logging.info("Filtering to unrelated individuals for PC calculations.")
            mt_norelateds = mt_filtered.filter_cols(mt_filtered.related_to_remove == False, keep=True)
            mt_norelateds = mt_norelateds.checkpoint(final_filtered, overwrite=True)

        else:
            logging.info("Detected final failing sample, variant, and genotype mt exists. Loading that.")
            mt_norelateds = hl.read_matrix_table(final_filtered)

        # MAF filter
        if (not utils.check_exists(final_maffilt)) or args.force:
            logging.info("Filtering to common variants and LD pruning dataset.")
            mt_mafpruned = vq.maf_filter(mt_norelateds, 0.05)
            mt_mafpruned = mt_mafpruned.checkpoint(final_maffilt)
        else:
            logging.info("Detected final MAF filtered mt exists. Loading that.")
            mt_mafpruned = hl.read_matrix_table(final_maffilt)

        # LD prune
        if (not utils.check_exists(final_ldpruned)) or args.force:
            logging.info('LD pruning final dataset for PC calculation')
            utils.remove_secondary(args.cluster_name, args.region)
            mt_ldpruned = vq.ld_prune(mt_mafpruned, args)
            mt_ldpruned = mt_ldpruned.checkpoint(final_ldpruned)
        else:
            logging.info("Detected final LDpruned mt exists. Loading that.")
            mt_ldpruned = hl.read_matrix_table(final_ldpruned)

        # Calculate PCs and project to relatives, plot
        mt = sq.project_pcs_relateds(mt_ldpruned, mt, args.pc_num)

        if args.pca_plot_annotations is not None:
            pca_annotations = args.pca_plot_annotations.strip().split(",")
            for annotation in pca_annotations:
                output_file(f"{datestr}_final_pcs_plot_{annotation}.html")
                p = hl.plot.scatter(mt.pc1, mt.pc2, label=mt[annotation])
                save(p)
        else:
            output_file(f'{datestr}_final_pcs_plot.html')
            pcplot = hl.plot.scatter(mt.pc1, mt.pc2)
            save(pcplot)

        # Export rows and columns
        mtcols = mt.cols()
        mtcols = mtcols.flatten()
        mtcols.export(os.path.join(args.out_dir, args.out_name + '_final_dataset_cols.tsv'))

        mtrows = mt.rows()
        mtrows = mtrows.flatten()
        mtrows.export(os.path.join(args.out_dir, args.out_name + '_final_dataset_rows.tsv'))

    # Send logs and finish-up notice
    logging.info('Pipeline ran successfully! Copying logs and shutting down cluster in 10 minutes.')
    utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

