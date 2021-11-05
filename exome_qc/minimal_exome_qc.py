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


    ##################################
    # Parse arguments for imput data #
    ##################################
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=str, required=True,
                        help="Name of VCF file (or files) to import, comma separated if > 1 file.")
    parser.add_argument('--cluster_name', type=str, help='Name of cluster for scaling in pipeline.')
    parser.add_argument("--region", type=str, default="europe-west1", help="Region of cluster for scaling.")
    parser.add_argument('--num_secondary_workers', type=int, default=20,
                        help='Number of secondary workers for scaling in applicable steps.')
    parser.add_argument("--out_file", type=str, help="Name of matrix table to output.")
    parser.add_argument("--log_dir", type=str, required=True, help="Location where logs should be written to.")
    parser.add_argument("--data_dir", type=str, required=True, help="Location where VCFs to import exist.")
    parser.add_argument("--out_dir", type=str, required=True, help="Location to write combined + vep annotated mt")
    parser.add_argument("--reference_genome", default='GRCh37', choices=['GRCh37', 'GRCh38'],
                        help="Reference genome build.")
    parser.add_argument("--chr_prefix", action='store_true', help="Chromosomes are of form 'chr1', NOT '1' etc.")
    parser.add_argument("--force_bgz", action='store_true', help="Force blog gzip import? Default true.")
    parser.add_argument("--call_fields", default="PGT", help="Name of genotype call field in VCF, default PGT.")
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

    args = parser.parse_args()

    hl.init()

    #####################################
    # Configure logging, define outputs #
    #####################################
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day
    log_file = 'import_minimal_qc-' + timestr + '.txt'

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

    ####################
    # Import VCF files #
    ####################
    vcf_files = args.vcf.strip().split(",")

    if (len(vcf_files) > 1) and (args.out_file is None):
        logging.error("Error! Must give matrix table file name with --out_file if importing more than one VCF.")
        exit(1)

    if args.test:
        test_str = "_test"
    else:
        test_str = ""

    if args.out_file is None:
        basename = os.path.basename(vcf_files[0]).replace(".vcf", "").replace(".gz", "").replace(".bgz", "")
    else:
        basename = args.out_file.rstrip("/").replace(".mt", "")

    out_basename = os.path.join(args.out_dir, basename)
    combined_mt_fn = out_basename + f"_combined{test_str}.mt/"
    split_fn = out_basename + f"_split{test_str}.mt/"

    #############################################
    # Combine VCF files and split multiallelics #
    #############################################
    if (not utils.check_exists(split_fn)) or args.force:

        # Import VCF files and combine
        if (not utils.check_exists(combined_mt_fn)) or args.force:
            mt = utils.load_vcfs(vcf_files, args.data_dir, args.out_dir, force=args.force, test=args.test,
                                 chr_prefix=args.chr_prefix, reference_genome=args.reference_genome,
                                 force_bgz=args.force_bgz,
                                 call_fields=args.call_fields)

            mt = mt.checkpoint(combined_mt_fn, overwrite=True)
            logging.info(f"Final matrix table count: {mt.count()}")
            utils.copy_logs_output(log_dir, log_file=log_file, plot_dir=args.data_dir)
        else:
            logging.info("Detected VCF file already converted to matrix table, skipping VCF import.")

        # Split multiallelics
        logging.info('Splitting multiallelic variants')
        mt = hl.read_matrix_table(combined_mt_fn)

        mt = hl.split_multi_hts(mt)
        mt = mt.checkpoint(split_fn, overwrite=True)
        logging.info('Split count: ' + str(mt.count()))
        utils.copy_logs_output(log_dir, log_file=log_file, plot_dir=args.data_dir)
    else:
        logging.info("Detected split mt exists, skipping splitting MT.")

    #######################
    # Low pass variant QC #
    #######################
    qcd_fn = out_basename + f"_low_pass_qcd{test_str}.mt/"

    if (not utils.check_exists(qcd_fn)) or args.force:
        logging.info("Running variant QC")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)
        mt = hl.read_matrix_table(split_fn)

        mt = vq.variant_quality_control(
            mt, qcd_fn, annotation_prefix="low_pass", min_dp=args.min_dp, min_gq=args.min_gq,
            max_het_ref_reads=args.max_het_ref_reads, min_het_ref_reads=args.min_het_ref_reads,
            min_hom_ref_ref_reads=args.min_hom_ref_ref_reads, max_hom_alt_ref_reads=args.max_hom_alt_ref_reads,
            call_rate=args.low_pass_min_call_rate, p_hwe=args.low_pass_p_hwe, snp_qd=args.snp_qd,
            indel_qd=args.indel_qd,
            ab_allowed_dev_het=args.ab_allowed_dev_het, count_failing=args.count_failing, sex_aware_call_rate=False,
            samples_qc=False, force=args.force
        )

        mt = mt.checkpoint(qcd_fn, overwrite=True)
    else:
        logging.info("Detected variant QC already performed, skipping that.")

    #############################
    # Export to VCF, text files #
    #############################
    # Filter out failing genotypes and variants, export to vcf
    mt_filt_fn = out_basename + f"_low_pass_filtered{test_str}.mt/"
    logging.info("Filtering out failing variants and genotypes and writing to VCF.")

    if (not utils.check_exists(mt_filt_fn)) or args.force:
        mt = hl.read_matrix_table(qcd_fn)

        mt_filt = sq.filter_failing(
            mt, mt_filt_fn, prefix='low_pass', entries=True, variants=True, samples=False, unfilter_entries=True,
            pheno_qc=False, min_dp=args.min_dp, min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
            min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
            max_hom_alt_ref_reads=args.max_hom_alt_ref_reads, force=args.force
        )

        mt_filt = mt_filt.checkpoint(mt_filt_fn, overwrite=True)

    else:
        mt_filt = hl.read_matrix_table(mt_filt_fn)
        mt = hl.read_matrix_table(qcd_fn)

    vcf_name = out_basename + "_failing_variants_genotypes_filtered.vcf.bgz"
    if (not utils.check_exists(vcf_name)) or args.force:
        hl.export_vcf(mt_filt, vcf_name, tabix=True)
    else:
        logging.info("Detected VCF already exported, skipping export.")

    # Export variant information to a separate tsv file
    variant_info_fn = out_basename + "_variant_info.tsv.bgz"

    if (not utils.check_exists(variant_info_fn)) or args.force:
        var_info = mt.rows()
        var_info = var_info.flatten()

        var_info.export(variant_info_fn)
    else:
        logging.info("Detected variant info table already exported, skipping export.")

