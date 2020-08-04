"""
Script for doing exome sequencing data quality control from a Hail matrix table input that has been VEP annotated.

Author: Lea Urpa, August 2020
"""
import argparse


def parse_arguments(arguments):
    """
    Takes arguments from the command line, runs argparse, and returns args object.

    :param arguments: command line arguments
    :return: args namespace object
    """
    parser = argparse.ArgumentParser(description="v9 exome sequencing dataset quality control pipeline.")

    # Pipeline parameters #
    params = parser.add_argument_group("Pipeline parameters")
    params.add_argument("--checkpoint", type=int, help="Checkpoint to start pipeline at.")
    params.add_argument("--test", action='store_true', help="run test with just chrom 22?")
    params.add_argument('--force', type=bool, default=True, help='Overwrite previous pipeline checkpoints?')
    params.add_argument('--skip_ab_filter', action='store_true', help='Skip allelic balance filtering?')
    params.add_argument('--run_king', action='store_true', help='Pause pipeline to run King relatedness calculations?')
    params.add_argument('--pc_num', default=10, help="Number of PCs to calculate.")
    params.add_argument('--verbosity', type=int, default=1,
                        help='Verbosity? Does counts, takes extra time for processing. 0 is least verbose.')
    params.add_argument('--num_preemptible_workers', type=int, default=100,
                        help='Number of preemptible workers for scaling in applicable steps.')

    # Pipeline inputs #
    inputs = parser.add_argument_group("Pipeline inputs and information.")
    inputs.add_argument("-file", type=str, help="Name of matrix table to run QC pipeline on.")
    inputs.add_argument("--out_name", type=str, help="Output name ")
    inputs.add_argument('--cluster_name', type=str, help='Name of cluster for scaling in pipeline.')
    inputs.add_argument("--out_dir", type=str, help="Directory to write output data to.")
    inputs.add_argument("--log_dir", type=str, help="Directory to write logs to.")
    inputs.add_argument("--scripts_dir", type=str, help="Directory containing python scripts.")
    inputs.add_argument("--bam_metadata", type=str, help="File containing bam metadata information.")
    inputs.add_argument("--samples_annotation_files", type=str,
                        help="Files to annotate the samples with, comma separated.")

    # Variant QC thresholds #
    var_thresh = parser.add_argument_group("Variant QC thresholds. If not indicated 'final' or 'low pass' in name, "
                                           "thresholds are used for both low pass variant filtering and final variant"
                                           "filtering.")
    var_thresh.add_argument("--low_pass_p_hwe", default=1e-9, help="Low pass variant QC HWE cutoff")
    var_thresh.add_argument("--final_p_hwe", default=1e-6, help="Final variant QC HWE cutoff")
    var_thresh.add_argument("--low_pass_min_call_rate", default=0.8, help="Low pass variant QC min call rate")
    var_thresh.add_argument("--final_min_call_rate", default=0.9, help="Final variant QC min call rate")
    var_thresh.add_argument("--min_dp", default=10, help="Variant QC min read depth")
    var_thresh.add_argument("--snp_qd", default=2, help="Variant QC min quality by depth for snps")
    var_thresh.add_argument("--indel_qd", default=3, help="Variant QC min quality by depth for indels")

    # Genotype QC thresholds #
    geno_thresh = parser.add_argument_group("Genotype QC thresholds. If not indicated 'final' or 'low pass' in name,"
                                            "thresholds are used for both low pass and final genotype filtering.")
    geno_thresh.add_argument("--min_het_ref_reads", default=0.2, help="min % reference reads for a het GT call")
    geno_thresh.add_argument("--max_het_ref_reads", default=0.8, help="max % reference reads for a het GT call")
    geno_thresh.add_argument("--min_hom_ref_ref_reads", default=0.9, help="min % reference reads for a ref GT call")
    geno_thresh.add_argument("--max_hom_alt_ref_reads", default=0.1, help="max % reference reads for an alt GT call")
    geno_thresh.add_argument("--low_pass_ab_allowed_dev_het", default=0.7,
                             help="% of het GT calls for a variant allowed to be out of allelic balance (% ref or alt "
                                  "reads out of range for het GT call), low pass genotype QC")
    geno_thresh.add_argument("--final_ab_allowed_dev_het", default=0.8,
                             help="% of het GT calls for a variant allowed to be out of allelic balance (% ref or alt "
                                  "reads out of range for het GT call), final genotype QC")

    # LD pruning thresholds #
    ld_thresh = parser.add_argument_group("Thresholds for LD pruning.")
    ld_thresh.add_argument("--r2", default=0.2, help="r2 correlation cutoff for LD pruning.")
    ld_thresh.add_argument("--bp_window_size", default=500000, help="Sliding window size for LD pruning in bp.")

    # Kinship thresholds #
    kin_thresh = parser.add_argument_group("Kinship thresholds.")
    kin_thresh.add_argument("--kinship_threshold", default=0.088, help="Kinship cutoff for finding related pairs.")
    kin_thresh.add_argument("--ind_maf", default=0.001, help="Minor allele frequency cutoff for calculating kinship.")
    kin_thresh.add_argument("--use_case_info", default=True,
                            help="Use case info to minimize case removal when finding unrelated set of individuals.")
    kin_thresh.add_argument("--plot_kin", default=True, help="Plot kinship values for visual inspection.")

    # Samples QC thresholds #
    samples_thresh = parser.add_argument_group("Samples QC thresholds.")
    samples_thresh.add_argument("--batches", action='store_true', help="Stratify samples QC by batch/cohort?")
    samples_thresh.add_argument("--batch_cohort_name", type=str, help="Samples annotation giving cohort or batch")
    samples_thresh.add_argument("--chimeras_max", default=0.05, help="Max % of chimeras allowed for sample")
    samples_thresh.add_argument("--contamination_max", default=0.05, help="Max % contamination allowed for sample")
    samples_thresh.add_argument("--std_dev", default=4,
                                help="Number of standard deviations from mean sample can deviate on heterozygosity.")
    samples_thresh.add_argument("--center_measure", default="mad", choices=['mad', 'std_dev'])

    # Impute sex thresholds #
    sex_thresh = parser.add_argument_group("Impute sex thresholds.")
    sex_thresh.add_argument("--female_threshold", default=0.4, help="F-stat cutoff for defining female")
    sex_thresh.add_argument("--male_threshold", default=0.8, help="F-stat cutoff for defining male")

    # Filter by phenotype thresholds #
    pheno_thresh = parser.add_argument_group("Thresholds for filtering by phenotype.")
    pheno_thresh.add_argument("--pheno_col", type=str, help="Samples annotation giving phenotype boolean annotation.")
    pheno_thresh.add_argument("--pheno_call_rate", default=0.95,
                              help="Min call rate for variant, in cases + controls separately.")

    parsed_args = parser.parse_args(arguments)

    return parsed_args


if __name__ == "__main__":
    ######################################
    # Initialize Hail and import scripts #
    ######################################
    import sys
    import os
    import logging
    import hail as hl
    hl.init()

    args = parse_arguments(sys.argv[1:])

    # Import python scripts to access helper functions #
    scripts = ["helper_scripts.py", "v9_exome_qc_parameters.py", "pipeline_functions.py"]
    for script in scripts:
        hl.spark_context().addPyFile(args.scripts_dir + script)

    import helper_scripts as h
    import pipeline_functions as p

    ####################
    # Configure logger #
    ####################
    logstem = 'import_vep_annotate-'
    datestr, timestr, log_file = h.configure_logging(logstem=logstem)

    log_dir = os.path.join(args.log_dir, logstem + datestr)

    # Configure logger
    root = logging.getLogger()
    log_formatter = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(filename=log_file, format=log_formatter, level=logging.INFO)

    handler = logging.StreamHandler(sys.stdout)
    root.addHandler(handler)

    ################
    # Run pipeline #
    ################
    args.checkpoint = 0
    args.output_stem = os.path.join(args.out_dir, args.out_name)
    args.checkpoint_folder = os.path.join(args.out_dir, "checkpoint_mts/")
    args.plot_folder = os.path.join(args.out_dir, "plots")

    try:
        # Load in data according to parameters given
        mt = p.load_data(args)

        # Annotate samples
        mt = p.annotate_samples(mt, args)

        # Low-pass variant QC
        mt = p.low_pass_var_qc(mt, args)

        # Phenotype Samples QC
        mt = p.phenotype_samples_qc(mt, args)

        # MAF pruning dataset
        mt, mt_mafpruned = p.maf_prune_relatedness(mt, args)

        # Export data to find related individuals in King, if necessary
        mt, mt_mafpruned = p.find_related_individuals(mt, mt_mafpruned, args)

        # LD prune MAF pruned dataset
        # This uses MAF 0.01 (1%) cutoff since that's what we exported to King before. Should be fine.
        mt, mt_ldpruned = p.ld_prune_popoutliers(mt, mt_mafpruned, args)

        # Find population outliers (excludes relatives)
        # (Excluding related individuals from analysis, but keeping them in the dataset)
        mt = p.find_pop_outliers(mt, mt_ldpruned, args)

        # Analytical samples QC (samples QC hard filters)
        # (Excluding population outliers from analysis, but keeping them in the dataset)
        mt = p.analytical_samples_qc(mt, args)

        # MAF prune dataset for sex imputation
        # This uses the default 0.05 (5%) MAF cutoff for common variants.
        mt, mt_mafpruned = p.maf_prune_sex_imputation(mt, args)

        # Impute sex
        # (Excluding analytical failing samples from analysis, but keeping them in the dataset)
        mt = p.impute_sex(mt, mt_mafpruned, args)

        # Variant QC filtering
        # (Excluding population outliers + analytical samples fails, but keeping them in the dataset)
        mt = p.variant_qc(mt, args)

        # Annotate gnomad + CADD
        mt = p.annotate_gnomad_cadd(mt, args)

        # Filter variants missing by pheno
        # (Excluding population outliers + analytical sample fails, but keeping them in the dataset)
        mt = p.filter_missing_by_pheno(mt, args)

        # Calculate final PCS
        # (Excluding population outliers + analytical sample fails + relatives but projecting them back into the PCS)
        mt = p.calculate_final_pcs(mt, args)

        # Get final case-control counts
        # (Excluding population outliers + analytical sample fails)
        mt = p.case_control_genotype_counts(mt,  args)

        # Export matrix table columns for optmatch
        mt = mt.filter_cols((mt.non_finns_to_remove == False) & (mt.fail_analytical == 0))
        mtcols = mt.cols()
        mtcols.export(args.output_stem + '_final_dataset_cols_passingQC.tsv')

        # Send logs and finish-up notice
        logging.info('Pipeline ran successfully! Copying logs and shutting down cluster in 10 minutes.')
        h.copy_logs_output(log_dir, timestr=timestr, log_file=log_file, plot_dir=args.plot_folder)

    except Exception as e:
        logging.error('Something went wrong! Copying logs and shutting down cluster in 10 minutes.')
        logging.error(e)
        h.copy_logs_output(log_dir, timestr=timestr, log_file=log_file, plot_dir=args.plot_folder)
