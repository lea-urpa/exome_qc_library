"""
Script for VCF import and minimal variant quality control.

Author: Lea Urpa, November 2021
"""

if __name__ == "__main__":
    import os
    import argparse
    import sys
    import logging
    import time
    import hail as hl
    import utils
    import variant_qc as vq
    import samples_qc as sq
    import samples_annotation as sa
    import variant_annotation as va

    ##################################
    # Parse arguments for imput data #
    ##################################
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=str, help="Name of VCF file (or files) to import, comma separated if > 1 file.")
    parser.add_argument("--mt", type=str, help="Name of (single) matrix table to load.")
    parser.add_argument("--vcfs_chrom_split", action="store_true", help="Are vcfs given chrom split, e.g. same samples?")
    parser.add_argument("--samples_annotation_files", type=str,
                        help="Files to annotate the samples with, comma separated.")
    parser.add_argument("--samples_col", type=str,
                        help="Name of samples column in sample annotation files. Must be the same in all files.")
    parser.add_argument("--samples_delim", type=str, default="\t",
                        help="Delimiter in sample annotation files. Must be the same in all files. Default tab.")
    parser.add_argument("--samples_miss", type=str, default="n/a",
                        help="String for missing values in annotation files, e.g. NA. Must be the same in all files.")
    parser.add_argument("--chimeras_col", type=str,
                                help="Column in matrix table or annotation files giving sample chimera percentage")
    parser.add_argument("--contamination_col", type=str,
                                help="Column in matrix table or annotation files giving sample contamination percent.")
    parser.add_argument("--batch_col_name", type=str,
                                help="Samples annotation in matrix table or annotation giving batch/cohort for "
                                     "stratified samples QC (TiTv, het/homvar, indel ratios, n singletons).")
    parser.add_argument("--sex_col", type=str, help="Column in annotation files givien sex column.")
    parser.add_argument("--male_str", type=str, help="string in sex column indicating a male sample.")
    parser.add_argument("--female_str", type=str, help="string in sex column indicating a female sample.")
    parser.add_argument('--cluster_name', type=str, help='Name of cluster for scaling in pipeline.')
    parser.add_argument("--region", type=str, default="europe-west1", help="Region of cluster for scaling.")
    parser.add_argument('--num_secondary_workers', type=int, default=20,
                        help='Number of secondary workers for scaling in applicable steps.')
    parser.add_argument("--out_file", type=str, help="Name of matrix table to output.")
    parser.add_argument("--log_dir", type=str, required=True, help="Location where logs should be written to.")
    parser.add_argument("--data_dir", type=str, help="Location where VCFs to import exist.")
    parser.add_argument("--out_dir", type=str, required=True, help="Location to write combined + vep annotated mt")
    parser.add_argument("--reference_genome", default='GRCh37', choices=['GRCh37', 'GRCh38'],
                        help="Reference genome build.")
    parser.add_argument("--split_by_chrom", action="store_true", help="Split vcf output by chromosome?")
    parser.add_argument("--chr_prefix", action='store_true', help="Chromosomes are of form 'chr1', NOT '1' etc.")
    parser.add_argument("--force_bgz", action='store_true', help="Force blog gzip import? Default true.")
    parser.add_argument("--call_fields", default="PGT", help="Name of genotype call field in VCF, default PGT.")
    parser.add_argument("--annotate_variants", action='store_true',
                        help="Annotate variants with LOF + missense genes they may lie in?")
    parser.add_argument("--do_not_export_vcf", action='store_false', help="Export QCd dataset as vcf file(s)?")
    parser.add_argument("--test", action='store_true', help="Filters data to just chr 22 for testing purposes.")
    parser.add_argument("--force", action='store_true', help="Force re-run of all steps?")

    # Genotype QC thresholds #
    geno_thresh = parser.add_argument_group("Genotype QC thresholds. If not indicated 'final' or 'low pass' in name,"
                                            "thresholds are used for both low pass and final genotype filtering.")
    geno_thresh.add_argument("--min_dp", default=10, type=float, help="min read depth for a genotype")
    geno_thresh.add_argument("--min_gq", default=20, type=float, help="min genotype quality for all GT calls.")
    geno_thresh.add_argument("--min_het_ref_reads", default=0.2, type=float,
                             help="min % reference reads for a het GT call")
    geno_thresh.add_argument("--max_het_ref_reads", default=0.8, type=float,
                             help="max % reference reads for a het GT call")
    geno_thresh.add_argument("--min_hom_ref_ref_reads", default=0.9, type=float,
                             help="min % reference reads for a ref GT call")
    geno_thresh.add_argument("--max_hom_alt_ref_reads", default=0.1, type=float,
                             help="max % reference reads for an alt GT call")
    geno_thresh.add_argument("--count_failing", default=True, type=bool,
                             help="Count number of genotypes failing each filter? Slow, but handy for troubleshooting.")

    # Variant QC thresholds #
    var_thresh = parser.add_argument_group("Variant QC thresholds.")
    var_thresh.add_argument("--low_pass_p_hwe", default=1e-9, type=float, help="Low pass variant QC HWE cutoff")
    var_thresh.add_argument("--low_pass_min_call_rate", default=0.8, type=float,
                            help="Low pass variant QC min call rate")
    var_thresh.add_argument("--snp_qd", default=2, type=float, help="Variant QC min quality by depth for snps")
    var_thresh.add_argument("--indel_qd", default=3, type=float, help="Variant QC min quality by depth for indels")
    var_thresh.add_argument("--ab_allowed_dev_het", default=0.8, type=float,
                            help="% of het GT calls for a variant that must be in allelic balance (% ref or alt "
                                 "reads out of range for het GT call)")

    # Samples QC thresholds #
    samples_thresh = parser.add_argument_group("Samples QC thresholds.")
    samples_thresh.add_argument("--skip_samples_qc", action="store_true", help="Skip samples QC?")
    samples_thresh.add_argument("--sample_call_rate", type=float, default=0.8,
                                help="Minimum genotype call rate per sample. Default none- samples not filtered on "
                                     "call rate.")
    samples_thresh.add_argument("--chimeras_max", default=0.05, type=float, help="Max % of chimeras allowed for sample")
    samples_thresh.add_argument("--contamination_max", default=0.05, type=float,
                                help="Max % contamination allowed for sample")
    samples_thresh.add_argument("--sampleqc_sd_threshold", default=4, type=int,
                                help="Number of standard deviations from mean sample can deviate on heterozygosity.")

    # Impute sex thresholds #
    sex_thresh = parser.add_argument_group("Impute sex thresholds.")
    sex_thresh.add_argument("--female_threshold", default=0.2, type=float, help="F-stat cutoff for defining female")
    sex_thresh.add_argument("--male_threshold", default=0.8, type=float, help="F-stat cutoff for defining male")

    args = parser.parse_args()

    hl.init()

    if not args.skip_samples_qc:
        if (args.chimeras_col is None) or (args.contamination_col is None):
            logging.error("Error! If not skipping samples QC, --chimeras_col and --contamination_col must be given.")
            exit(1)

    #####################################
    # Configure logging, define outputs #
    #####################################
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day
    log_file = 'import_minimal_qc-' + timestr + '.txt'
    plot_dir = args.log_dir.rstrip("/") + "/plots"

    root = logging.getLogger()  # creates logger
    root.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Add file handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add streaming handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    log_dir = os.path.join(args.log_dir, 'import_vep_annotate-' + datestr)

    #########################################
    # Check dataproc cluster vs file inputs #
    #########################################
    for file_url in [args.log_dir, args.out_dir, ]:
        utils.check_regions(args.region, file_url)

    if args.test:
        test_str = "_test"
    else:
        test_str = ""

    #######################################
    # Import VCF files, combine and split #
    #######################################
    if args.vcf is not None:
        vcf_files = args.vcf.strip().split(",")

        if (len(vcf_files) > 1) and (args.out_file is None):
            logging.error("Error! Must give matrix table file name with --out_file if importing more than one VCF.")
            exit(1)

        if args.data_dir is None:
            logging.error("Error! Directory of VCF files must be given with --data_dir.")
            exit(1)
        else:
            utils.check_regions(args.region, args.data_dir)

        if args.out_file is None:
            basename = os.path.basename(vcf_files[0]).replace(".vcf", "").replace(".gz", "").replace(".bgz", "")
        else:
            basename = args.out_file.rstrip("/").replace(".mt", "")

        out_basename = os.path.join(args.out_dir, basename)

        combined_mt_fn = out_basename + f"_combined{test_str}.mt/"
        split_fn = out_basename + f"_split{test_str}.mt/"

        if (not utils.check_exists(split_fn)) or args.force:
            # Import VCF files and combine
            if (not utils.check_exists(combined_mt_fn)) or args.force:
                mt = utils.load_vcfs(vcf_files, args.data_dir, args.out_dir, force=args.force, test=args.test,
                                     chr_prefix=args.chr_prefix, reference_genome=args.reference_genome,
                                     force_bgz=args.force_bgz,
                                     call_fields=args.call_fields, chrom_split=args.vcfs_chrom_split)

                mt = mt.checkpoint(combined_mt_fn, overwrite=True)
                logging.info(f"Final matrix table count: {mt.count()}")
                utils.copy_logs_output(log_dir, log_file=log_file, plot_dir=plot_dir)
            else:
                logging.info("Detected VCF file already converted to matrix table, skipping VCF import.")

            # Split multiallelics
            logging.info('Splitting multiallelic variants')
            mt = hl.read_matrix_table(combined_mt_fn)

            mt = hl.split_multi_hts(mt)
            mt = mt.checkpoint(split_fn, overwrite=True)
            logging.info('Split count: ' + str(mt.count()))
            utils.copy_logs_output(log_dir, log_file=log_file, plot_dir=plot_dir)
        else:
            logging.info("Detected split mt exists, skipping splitting MT.")

    elif args.mt is not None:
        utils.check_regions(args.region, args.mt)
        split_fn = args.mt
        basename = args.out_file.rstrip("/").replace(".mt", "")
        out_basename = os.path.join(args.out_dir, basename)
        if not utils.check_exists(args.mt):
            logging.error(f"Error! file {args.mt} does not exist! Exiting now.")
        else:
            logging.info(f"Loading mt {args.mt} for quality control.")

    else:
        logging.error("Error! Either matrix table or vcf must be given as input. Exiting now.")
        exit(0)

    #######################
    # Low pass variant QC #
    #######################
    counter = 1
    variant_qcd_fn = f"{args.out_dir}{counter}_{basename}_low_pass_qcd{test_str}.mt/"

    if (not utils.check_exists(variant_qcd_fn)) or args.force:
        logging.info("Running variant QC")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)
        mt = hl.read_matrix_table(split_fn)

        mt = vq.variant_quality_control(
            mt, variant_qcd_fn, annotation_prefix="low_pass", min_dp=args.min_dp, min_gq=args.min_gq,
            max_het_ref_reads=args.max_het_ref_reads, min_het_ref_reads=args.min_het_ref_reads,
            min_hom_ref_ref_reads=args.min_hom_ref_ref_reads, max_hom_alt_ref_reads=args.max_hom_alt_ref_reads,
            call_rate=args.low_pass_min_call_rate, p_hwe=args.low_pass_p_hwe, snp_qd=args.snp_qd,
            indel_qd=args.indel_qd,
            ab_allowed_dev_het=args.ab_allowed_dev_het, count_failing=args.count_failing, sex_aware_call_rate=False,
            samples_qc=False, force=args.force
        )

        mt = mt.checkpoint(variant_qcd_fn, overwrite=True)
        args.force = True
    else:
        logging.info("Detected variant QC already performed, skipping that.")

    ##############
    # Impute sex #
    ##############
    counter += 1
    sex_imputed = f"{args.out_dir}{counter}_{basename}_sex_imputed{test_str}.mt"
    filtered_nohwe = out_basename + f"_variant_filtered_nohwe{test_str}_tmp.mt"
    filtered_annot = out_basename + f"_variant_filtere_sex_annotated{test_str}_tmp.mt"

    if (not utils.check_exists(sex_imputed)) or args.force:
        if args.annotate_variants:
            logging.info("Annotating variant types:")

            mt = va.annotate_variants(mt)

        logging.info("Imputing sex")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        mt = hl.read_matrix_table(variant_qcd_fn)

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

            mt_count = mt_filtered.count_rows()
            if not mt_count > 20000:
                logging.info(f"Error! Not enough variants after filtering to passing all QC measures. "
                             f"var count: {mt_count}")
                exit()
            mt_filtered = mt_filtered.checkpoint(filtered_nohwe, overwrite=True)
            args.force = True
        else:
            mt_filtered = hl.read_matrix_table(filtered_nohwe)

        # Impute sex, annotate to main matrix table and filtered matrix table
        imputed_sex = sq.impute_sex_plot(mt_filtered, female_threshold=args.female_threshold,
                                         male_threshold=args.male_threshold, aaf_threshold=0.05)

        mt_filtered = mt_filtered.annotate_cols(is_female_imputed=imputed_sex[mt_filtered.s].is_female)
        mt = mt.annotate_cols(is_female_imputed=imputed_sex[mt.s].is_female, f_stat=imputed_sex[mt.s].f_stat)
        mt = mt.annotate_globals(
            sex_imputation_thresholds={'female_threshold': args.female_threshold,
                                       'male_threshold': args.male_threshold})

        # Annotate sex-aware variant annotations (gt filt only to have for all variants)
        gt_filt_fn = sex_imputed.rstrip("/").replace(".mt", "") + "_GT_filtered.mt/"
        mt_gt_filt = hl.read_matrix_table(gt_filt_fn)
        mt_gt_filt = mt_gt_filt.annotate_cols(is_female_imputed=imputed_sex[mt_gt_filt.s].is_female)
        mt_gt_filt, annotations_to_transfer = va.sex_aware_variant_annotations(mt_gt_filt)

        # Annotate sex-aware sample annotations, checkpoint
        mt_filtered = sa.sex_aware_sample_annotations(mt_filtered)
        for annotation in annotations_to_transfer:
            mt_filtered = mt_filtered.annotate_rows(**{annotation: mt_gt_filt.rows()[mt_filtered.row_key][annotation]})
        mt_filtered = mt_filtered.checkpoint(filtered_annot, overwrite=True)

        # Annotate main MT with variant and sample sex-aware annotations
        for annotation in annotations_to_transfer:
            mt = mt.annotate_rows(**{annotation: mt_filtered.rows()[mt.row_key][annotation]})
        mt = mt.annotate_cols(sexaware_sample_call_rate=mt_filtered.cols()[mt.s].sexaware_sample_call_rate)

        mt = mt.checkpoint(sex_imputed, overwrite=True)
    else:
        logging.info("Detected sex imputed, skipping sex imputation.")

    ##################
    # Run samples QC #
    ##################
    if not args.skip_samples_qc:
        counter += 1
        samples_qcd_fn = f"{args.out_dir}{counter}_{basename}_samples_qcd{test_str}.mt/"
        annotated_fn = out_basename + f"_samples_qcd_samples_annotated{test_str}_tmp.mt"
        filtered_mt_fn = out_basename + f"_samples_qcd_gts_vars_filtered{test_str}_tmp.mt"
        first_samples_qc_fn = out_basename + f"_samples_qcd_no_sex_check{test_str}_tmp.mt"

        if (not utils.check_exists(samples_qcd_fn)) or args.force:
            mt = hl.read_matrix_table(sex_imputed)

            # Annotate with bam metadata
            if (not utils.check_exists(annotated_fn)) or args.force:
                mt = sa.annotate_cols_from_file(mt, args.samples_annotation_files, args.samples_delim, args.samples_col,
                                                args.samples_miss)
                mt = mt.checkpoint(annotated_fn, overwrite=True)
                args.force = True
            else:
                mt = hl.read_matrix_table(annotated_fn)

            # Check columns exist
            for colname in ['chimeras_col', 'contamination_col', 'sex_col']:
                col = getattr(args, colname)
                try:
                    test = hl.is_defined(mt[col])
                except Exception as e:
                    logging.error(f"Error! Given column annotation {col} does not actually exist after inputting sample "
                                  f"annotations.")
                    logging.error(e)
                    exit(1)

            # Filter failing variants and genotypes
            if (not utils.check_exists(filtered_mt_fn)) or args.force:
                logging.info("Filtering out failing genotypes and variants.")
                mt_filtered = sq.filter_failing(
                    mt, samples_qcd_fn, prefix='low_pass', entries=True, variants=True, samples=False, unfilter_entries=False,
                    pheno_qc=False, min_dp=args.min_dp, min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
                    min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
                    max_hom_alt_ref_reads=args.max_hom_alt_ref_reads, force=args.force
                )
                mt_filtered = mt_filtered.checkpoint(filtered_mt_fn, overwrite=True)
                args.force = True
            else:
                mt_filtered = hl.read_matrix_table(filtered_mt_fn)

            # Run samples QC
            if (not utils.check_exists(first_samples_qc_fn)) or args.force:
                mt = sq.samples_qc(
                    mt_filtered, mt, samples_qcd_fn, count_failing=args.count_failing, sample_call_rate=args.sample_call_rate,
                    chimeras_col=args.chimeras_col, chimeras_max=args.chimeras_max, contamination_col=args.contamination_col,
                    contamination_max=args.contamination_max, batch_col_name=args.batch_col_name,
                    sampleqc_sd_threshold=args.sampleqc_sd_threshold, force=args.force
                )
                mt = mt.checkpoint(first_samples_qc_fn, overwrite=True)
                args.force = True
            else:
                mt = hl.read_matrix_table(first_samples_qc_fn)

            # Convert sex_col to is_female column
            mt = mt.annotate_cols(is_female_reported=hl.cond(
                hl.is_defined(mt[args.sex_col]) & (mt[args.sex_col] == args.male_str),
                False,
                hl.cond(
                    hl.is_defined(mt[args.sex_col]) & (mt[args.sex_col] == args.female_str),
                    True, hl.null(hl.tbool)
                )
            ))

            # Check that imputed sex matches given sex
            mt = mt.annotate_cols(failing_samples_qc=hl.cond(
                hl.is_defined(mt.is_female_imputed) & hl.is_defined(mt[args.sex_col]) &
                (mt.is_female_imputed != mt.is_female_reported),
                mt.failing_samples_qc.append("missing_sexaware_sample_call_rate"),
                mt.failing_samples_qc
            ))

            logging.info(f"Writing checkpoint: sample QC")
            mt = mt.checkpoint(samples_qcd_fn, overwrite=True)
            utils.copy_logs_output(args.log_dir, log_file=log_file, plot_dir=plot_dir)
        else:
            logging.info("Detected samples QC completed, skipping this step.")
    else:
        logging.info("User indicated samples QC should be skipped. Moving on.")
        samples_qcd_fn = sex_imputed

    # Add pop oulier analysis

    #####################################
    # Filter failing gts, vars, samples #
    #####################################
    # Filter out failing genotypes, samples, and variants, export to vcf
    counter += 1
    mt_filt_fn = f"{args.out_dir}{counter}_{basename}_failing_filtered{test_str}.mt/"
    logging.info("Filtering out failing variants, genotypes, and samples to write to VCF.")

    if (not utils.check_exists(mt_filt_fn)) or args.force:
        mt = hl.read_matrix_table(samples_qcd_fn)

        if args.skip_samples_qc:
            filter_samples = False
        else:
            filter_samples = True
        mt_filt = sq.filter_failing(
            mt, mt_filt_fn, prefix='low_pass', entries=True, variants=True, samples=filter_samples, unfilter_entries=True,
            pheno_qc=False, min_dp=args.min_dp, min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
            min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
            max_hom_alt_ref_reads=args.max_hom_alt_ref_reads, force=args.force, pop_outliers=False
        )

        # Add final variant QC measures after filtering
        mt_filt = hl.variant_qc(mt_filt, name='final_variant_qc_')

        mt_filt = mt_filt.checkpoint(mt_filt_fn, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=log_file, plot_dir=plot_dir)
        args.force = True
    else:
        mt_filt = hl.read_matrix_table(mt_filt_fn)
        mt = hl.read_matrix_table(samples_qcd_fn)

    ################################################################
    # Export unfiltered variant information to a separate tsv file #
    ################################################################
    logging.info("Exporting sample (column) and variant (row) data.")
    variant_info_fn = out_basename + "_unfiltered_variant_info.tsv.bgz"

    if (not utils.check_exists(variant_info_fn)) or args.force:
        args.force = True
        var_info = mt.rows()
        var_info = var_info.flatten()

        var_info.export(variant_info_fn)
    else:
        logging.info("Detected variant info table already exported, skipping export.")

    # Export samples information to a separate tsv file
    sample_info_fn = out_basename + "_unfiltered_sample_info.tsv.bgz"

    if (not utils.check_exists(sample_info_fn)) or args.force:
        args.force = True
        sample_info = mt.cols()
        sample_info = sample_info.flatten()

        sample_info.export(sample_info_fn)
    else:
        logging.info("Detected sample info table already exported, skipping export.")

    utils.copy_logs_output(args.log_dir, log_file=log_file, plot_dir=plot_dir)

    #############################################
    # Export filtered mt as VCF, per chromosome #
    #############################################
    if not args.do_not_export_vcf:
        utils.remove_secondary(args.cluster_name, region=args.region)
        if args.split_by_chrom:
            logging.info("Exporting VCFs split by chromosome")
            chroms = [str(x) for x in range(1, 23)]
            chroms.extend(["X", "Y"])
            if args.reference_genome == "GRCh38":
                chroms = [f"chr{x}" for x in chroms]

            for chrom in chroms:
                if args.test:
                    if not "22" in chrom:
                        pass

                vcf_name = out_basename + f"_failing_filtered_chrom_{chrom}.vcf.bgz"
                if (not utils.check_exists(vcf_name)) or args.force:

                    mt_tmp = mt_filt.filter_rows(mt_filt.locus.contig == chrom)
                    if mt_tmp.count_rows() > 0:
                        hl.export_vcf(mt_tmp, vcf_name, tabix=True)
                        args.force = True
                    else:
                        logging.warning(f"Warning! No variants for chromosome {chrom} in dataset!")
                else:
                    logging.info(f"Detected {os.path.basename(vcf_name)} already exported, skipping export.")
            utils.copy_logs_output(args.log_dir, log_file=log_file, plot_dir=plot_dir)
        else:
            logging.info("Exporting VCF as one file.")
            vcf_name = out_basename + f"_failing_filtered.vcf.bgz"
            if (not utils.check_exists(vcf_name)) or args.force:
                args.force = True
                hl.export_vcf(mt_filt, vcf_name, tabix=True)
            else:
                logging.info(f"Detected {os.path.basename(vcf_name)} already exported, skipping export.")

            utils.copy_logs_output(args.log_dir, log_file=log_file, plot_dir=plot_dir)


