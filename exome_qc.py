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
    params.add_argument("--reference_genome", type=str, help="Reference_genome", choices=["GRCh37", "GRCh38"],
                        default="GRCh38")
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
    inputs.add_argument("-mt", type=str, help="Name of matrix table to run QC pipeline on.")
    inputs.add_argument("--out_name", type=str, help="Output name ")
    inputs.add_argument('--cluster_name', type=str, help='Name of cluster for scaling in pipeline.')
    inputs.add_argument("--out_dir", type=str, help="Directory to write output data to.")
    inputs.add_argument("--log_dir", type=str, help="Directory to write logs to.")
    inputs.add_argument("--scripts_dir", type=str, help="Directory containing python scripts.")
    inputs.add_argument("--bam_metadata", type=str, help="File containing bam metadata information.")
    inputs.add_argument("--samples_annotation_files", type=str,
                        help="Files to annotate the samples with, comma separated.")
    inputs.add_argument("--samples_col", type=str,
                        help="Name of samples column in sample annotation files. Must be the same in all files.")
    inputs.add_argument("--samples_delim", type=str,
                        help="Delimiter in sample annotation files. Must be the same in all files.")
    inputs.add_argument("--samples_miss", type=str,
                        help="String for missing values in annotation files, e.g. NA. Must be the same in all files.")
    inputs.add_argument("--fam_id", type=str, help="column name corresponding to sample's family ID. Used in kinship.")
    inputs.add_argument("--pat_id", type=str, help="column name corresponding to sample's paternal ID. Used in kinship.")
    inputs.add_argument("--mat_id", type=str, help="column name corresponding to sample's maternal ID. Used in kinship.")

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
    geno_thresh.add_argument("--min_gq", default=20, help="min genotype quality for all GT calls.")
    geno_thresh.add_argument("--min_het_ref_reads", default=0.2, help="min % reference reads for a het GT call")
    geno_thresh.add_argument("--max_het_ref_reads", default=0.8, help="max % reference reads for a het GT call")
    geno_thresh.add_argument("--min_hom_ref_ref_reads", default=0.9, help="min % reference reads for a ref GT call")
    geno_thresh.add_argument("--max_hom_alt_ref_reads", default=0.1, help="max % reference reads for an alt GT call")
    geno_thresh.add_argument("--ab_allowed_dev_het", default=0.8,
                             help="% of het GT calls for a variant allowed to be out of allelic balance (% ref or alt "
                                  "reads out of range for het GT call)")

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
    kin_thresh.add_argument("--relatives_removal_file", type=str,
                            help="File of related individuals to remove, one sample per line.")
    #TODO run check at beginning to see that relatives file is given if run_king is False

    # Samples removal options #
    samples_removal = parser.add_argument_group("Arbitrary sample removal options.")
    samples_removal.add_argument("--sample_removal_strings", type=str,
                                 help="string or strings that indicate a sample should be filtered out of the dataset."
                                      "If more than one, comma separated.")
    samples_removal.add_argument("--sample_removal_list", type=str,
                                 help="File containing one sample ID per line that should be removed from the dataset.")

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
    pheno_thresh.add_argument("--pheno_col", type=str,
                              help="Samples annotation giving phenotype boolean annotation. Note: samples missing"
                                   "phenotype information are ignored in many cases.")
    pheno_thresh.add_argument("--pheno_call_rate", default=0.95,
                              help="Min call rate for variant, in cases + controls separately.")

    # TODO add check to make sure all input files exist before running pipeline

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
    import qc_pipeline_functions as qc

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

    ##################################################
    # Pre-define checkpoint names to load data later #
    ##################################################
    args.sample_annotation = {'prefix': '1', 'suffix': '_after_sample_annotations'}
    args.low_pass_variant_qc = {'prefix': '2', 'suffix': '_after_low_pass_variant_qc'}

    ################
    # Run pipeline #
    ################
    args.cpcounter = 0
    args.lowpass_fail_name = 'failing_lowpass_varqc'
    args.final_fail_name = 'failing_final_varqc'
    args.output_stem = os.path.join(args.out_dir, args.out_name)
    args.checkpoint_folder = os.path.join(args.out_dir, "checkpoint_mts/")
    args.plot_folder = os.path.join(args.out_dir, "plots")

    try:
        # Load in data according to parameters given
        mt = qc.load_data(args)

        # Annotate samples
        mt = qc.annotate_samples(mt, args)

        # Phenotype Samples QC
        mt = qc.remove_samples(mt, args)

        # Low-pass variant QC
        mt = qc.low_pass_var_qc(mt, args)

        # Export data to find related individuals in King, if necessary
        mt, mt_mafpruned = qc.find_related_individuals(mt, args)

        # LD prune MAF pruned dataset
        # This uses MAF 0.01 (1%) cutoff since that's what we exported to King before. Should be fine.
        mt, mt_ldpruned = qc.ld_prune_popoutliers(mt, mt_mafpruned, args)

        # Find population outliers (excludes relatives)
        # (Excluding related individuals from analysis, but keeping them in the dataset)
        mt = qc.find_pop_outliers(mt, mt_ldpruned, args)

        # Analytical samples QC (samples QC hard filters)
        # (Excluding population outliers from analysis, but keeping them in the dataset)
        mt = qc.samples_qc(mt, args)

        # MAF prune dataset for sex imputation
        # This uses the default 0.05 (5%) MAF cutoff for common variants.
        mt, mt_mafpruned = qc.maf_prune_sex_imputation(mt, args)

        # Impute sex
        # (Excluding analytical failing samples from analysis, but keeping them in the dataset)
        mt = qc.impute_sex(mt, mt_mafpruned, args)

        # TODO make sure sex_aware_variant_annotations run before this step, or within this step before var qc
        # Variant QC filtering
        # (Excluding population outliers + analytical samples fails, but keeping them in the dataset)
        mt = qc.variant_qc(mt, args)

        # Annotate gnomad + CADD
        mt = qc.annotate_gnomad_cadd(mt, args)

        # Filter variants missing by pheno
        # (Excluding population outliers + analytical sample fails, but keeping them in the dataset)
        mt = qc.filter_missing_by_pheno(mt, args)

        # Calculate final PCS
        # (Excluding population outliers + analytical sample fails + relatives but projecting them back into the PCS)
        mt = qc.calculate_final_pcs(mt, args)

        # Get final case-control counts
        # (Excluding population outliers + analytical sample fails)
        mt = qc.case_control_genotype_counts(mt,  args)

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
