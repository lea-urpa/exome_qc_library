"""
Standalone pipeline to find duplicates in an exome sequencing dataset.

Author: Lea Urpa, October 2021
"""

import sys
import os
import logging
import time
import argparse
import hail as hl
from bokeh.io import output_file, save
import utils
import variant_qc as vq
import samples_qc as sq


if __name__ == "__main__":

    hl.init()

    ###################
    # Parse arguments #
    ###################
    # TODO add missing args from pipe below
    parser = argparse.ArgumentParser(description="Pipeline to find duplicates from an exome sequencing dataset.")
    parser.add_argument("--vcf", type=str, required=True,
                        help="Name of VCF file (or files) to import, comma separated if > 1 file.")
    parser.add_argument("--out_file", type=str, help="Name of matrix table to output.")
    parser.add_argument("--log_dir", type=str, required=True, help="Location where logs should be written to.")
    parser.add_argument("--data_dir", type=str, required=True, help="Location where VCFs to import exist.")
    parser.add_argument("--out_dir", type=str, required=True, help="Location to write combined + vep annotated mt")
    parser.add_argument("--reference_genome", default='GRCh38', choices=['GRCh37', 'GRCh38'],
                        help="Reference genome build.")
    parser.add_argument("--chr_prefix", action='store_true', help="Chromosomes are of form 'chr1', NOT '1' etc.")
    parser.add_argument("--force_bgz", action='store_true', help="Force blog gzip import? Default true.")
    parser.add_argument("--call_fields", default="PGT", help="Name of genotype call field in VCF, default PGT.")
    parser.add_argument("--test", action='store_true', help="Filters data to just chr 22 for testing purposes.")
    parser.add_argument("--force", action='store_true', help="Force re-run of all steps?")
    parser.add_argument('--cluster_name', type=str, help='Name of cluster for scaling in pipeline.')
    parser.add_argument('--num_secondary_workers', type=int, default=20,
                        help='Number of secondary workers for scaling in applicable steps.')
    parser.add_argument("--region", default='europe-west1', help='Region name for scaling in pipeline.')

    # Variant QC thresholds #
    var_thresh = parser.add_argument_group("Variant QC thresholds. If not indicated 'final' or 'low pass' in name, "
                                           "thresholds are used for both low pass variant filtering and final variant"
                                           "filtering.")
    var_thresh.add_argument("--low_pass_p_hwe", default=1e-9, type=float, help="Low pass variant QC HWE cutoff")
    var_thresh.add_argument("--low_pass_min_call_rate", default=0.8, type=float,
                            help="Low pass variant QC min call rate")
    var_thresh.add_argument("--snp_qd", default=2, type=float, help="Variant QC min quality by depth for snps")
    var_thresh.add_argument("--indel_qd", default=3, type=float, help="Variant QC min quality by depth for indels")
    var_thresh.add_argument("--ab_allowed_dev_het", default=0.8, type=float,
                            help="% of het GT calls for a variant that must be in allelic balance (% ref or alt "
                                 "reads out of range for het GT call)")

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

    args = parser.parse_args()

    # LD pruning thresholds #
    ld_thresh = parser.add_argument_group("Thresholds for LD pruning.")
    ld_thresh.add_argument("--r2", default=0.2, type=float, help="r2 correlation cutoff for LD pruning.")
    ld_thresh.add_argument("--bp_window_size", default=500000, type=int,
                           help="Sliding window size for LD pruning in bp.")

    # Kinship thresholds #
    kin_thresh = parser.add_argument_group("Kinship thresholds.")
    kin_thresh.add_argument("--ind_maf", default=0.05, type=float,
                            help="Minor allele frequency cutoff for calculating kinship.")
    kin_thresh.add_argument("--kinship_threshold", default=0.0883,
                            help="Threshold for kinship coefficient, above which individuals are defined as related.")

    ###############################################
    # Import scripts, configure logger and inputs #
    ###############################################
    hl.init()

    ## Configure logger ##
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day
    args.log_file = 'find_duplicates_' + timestr + '.txt'

    root = logging.getLogger()  # creates logger
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
    if args.test:
        args.test_str = "_test"
    else:
        args.test_str = ""

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

    if (not utils.check_exists(combined_mt_fn)) or args.force:
        mt = utils.load_vcfs(vcf_files, args.data_dir, force=args.force, test=args.test,
                             chr_prefix=args.chr_prefix, reference_genome=args.reference_genome,
                             force_bgz=args.force_bgz,
                             call_fields=args.call_fields)

        mt = mt.checkpoint(combined_mt_fn, overwrite=True)
    else:
        mt = hl.read_matrix_table(combined_mt_fn)

    logging.info(f"Final matrix table count: {mt.count()}")

    #######################
    # Split multiallelics #
    #######################
    split_fn = out_basename + f"_split{test_str}.mt/"

    if (not utils.check_exists(split_fn)) or args.force:
        logging.info('Splitting multiallelic variants')
        mt_split = hl.split_multi_hts(mt)
        mt_split = mt_split.checkpoint(split_fn, overwrite=True)
    else:
        logging.info("Detected split mt exists, loading that.")
        mt_split = hl.read_matrix_table(split_fn)
    logging.info('Split count: ' + str(mt_split.count()))

    #######################
    # low-pass variant QC #
    #######################
    low_pass_qcd = os.path.join(args.out_dir, f"{args.out_name}_low_pass_qcd{args.test_str}.mt/")

    if (not utils.check_exists(low_pass_qcd)) or args.force:
        logging.info("Running low-pass variant QC and genotype QC")

        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)
        mt = hl.read_matrix_table(split_fn)

        mt = vq.variant_quality_control(
            mt, low_pass_qcd, annotation_prefix="low_pass", min_dp=args.min_dp, min_gq=args.min_gq,
            max_het_ref_reads=args.max_het_ref_reads, min_het_ref_reads=args.min_het_ref_reads,
            min_hom_ref_ref_reads=args.min_hom_ref_ref_reads, max_hom_alt_ref_reads=args.max_hom_alt_ref_reads,
            call_rate=args.low_pass_min_call_rate, p_hwe=args.low_pass_p_hwe, snp_qd=args.snp_qd,
            indel_qd=args.indel_qd,
            ab_allowed_dev_het=args.ab_allowed_dev_het, count_failing=args.count_failing, sex_aware_call_rate=False,
            pheno_col=args.pheno_col, samples_qc=False, force=args.force
        )

        logging.info(f"Writing checkpoint after low pass variant QC")
        mt = mt.checkpoint(low_pass_qcd, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
    else:
        logging.info("Detected low-pass variant QC mt exists, skipping low-pass variant QC.")

    #########################
    # Calculate relatedness #
    #########################
    relatedness_calculated = os.path.join(args.out_dir, f"{args.out_name}_relatedness_calculated{args.test_str}.mt/")

    ld_pruned = os.path.join(args.out_dir, f"{args.out_name}_ld_pruned{args.test_str}.mt/")
    ld_pruned_maffilt = os.path.join(args.out_dir, f"{args.out_name}_maf_filt{args.test_str}.mt/")
    ld_pruned_annot = os.path.join(args.out_dir, f"{args.out_name}_ld_pruned_related{args.test_str}.mt/")

    if (not utils.check_exists(relatedness_calculated)) or args.force:
        logging.info("Calculating relatedness")
        mt = hl.read_matrix_table(low_pass_qcd)
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        ## LD prune and checkpoint ##
        if (not utils.check_exists(ld_pruned)) or args.force:

            # Filter failing samples, variants, and genotypes
            mt_gt_filt = sq.filter_failing(
                mt, ld_pruned, prefix='low_pass', entries=True, variants=False, samples=False, unfilter_entries=True,
                pheno_qc=False, min_dp=args.min_dp, min_gq=args.min_gq, max_het_ref_reads=args.max_het_ref_reads,
                min_het_ref_reads=args.min_het_ref_reads, min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
                max_hom_alt_ref_reads=args.max_hom_alt_ref_reads, force=args.force
            )

            mt_filtered = mt_gt_filt.filter_rows(
                mt_gt_filt.low_pass_failing_variant_qc.contains("failing_QD") |
                mt_gt_filt.low_pass_failing_variant_qc.contains("failing_VQSR_filters") |
                mt_gt_filt.low_pass_failing_variant_qc.contains("failing_call_rate"), keep=False
            )

            # Filter out low MAF variants
            if (not utils.check_exists(ld_pruned_maffilt)) or args.force:
                mt_maffilt = vq.maf_filter(mt_filtered, args.ind_maf)
                mt_maffilt = mt_maffilt.checkpoint(ld_pruned_maffilt, overwrite=True)
            else:
                mt_maffilt = hl.read_matrix_table(ld_pruned_maffilt)

            # LD prune if row count >80k
            mt_ldpruned = vq.downsample_variants(
                mt_maffilt, 80000, ld_pruned, r2=args.r2, bp_window_size=args.bp_window_size, ld_prune=True)

            logging.info(f"Writing checkpoint after LD pruned dataset")
            mt_ldpruned = mt_ldpruned.checkpoint(ld_pruned, overwrite=True)
        else:
            logging.info("Detected LD pruned dataset written, loading that.")
            mt_ldpruned = hl.read_matrix_table(ld_pruned)

        ## Calculate relatedness with King ##
        # THIS triggers shuffles, think about removing secondaries... seemed to run ok though
        if args.reference_genome is "GRCh38":
            autosomes = ["chr" + str(i) for i in range(1, 23)]
        else:
            autosomes = [str(i) for i in range(1, 23)]

        mt_autosomes = mt_ldpruned.filter_rows(hl.literal(autosomes).contains(mt_ldpruned.locus.contig))

        related_to_remove, related_info_ht = sq.king_relatedness(
            mt_autosomes, relatedness_calculated, kinship_threshold=args.kinship_threshold, pheno_col=args.pheno_col,
            force=args.force, cluster_name=args.cluster_name, num_secondary_workers=args.num_secondary_workers,
            region=args.region, reference_genome=args.reference_genome, export_duplicates=True)

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

        mt_ldpruned = mt_ldpruned.annotate_cols(
            related_num_connections=hl.or_else(mt_ldpruned.related_num_connections, 0))

        logging.info(f"Writing checkpoint after relatedness annotated")
        mt = mt.checkpoint(relatedness_calculated, overwrite=True)
        mt_ldpruned = mt_ldpruned.checkpoint(ld_pruned_annot, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)

    else:
        logging.info("Detected mt with relatives annotated exists, skipping relatedness calculation.")

