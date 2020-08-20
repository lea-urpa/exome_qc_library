"""
Script for doing exome sequencing data quality control from a Hail matrix table input that has been VEP annotated.

Author: Lea Urpa, August 2020
"""
import subprocess
import shlex
import argparse
import time


def parse_arguments(arguments):
    """
    Takes arguments from the command line, runs argparse, and returns args object.

    :param arguments: command line arguments
    :return: args namespace object
    """
    parser = argparse.ArgumentParser(description="v9 exome sequencing dataset quality control pipeline.")

    # Pipeline parameters #
    params = parser.add_argument_group("Pipeline parameters")
    params.add_argument("--checkpoint", type=int, help="Checkpoint to start pipeline at, default start (0).", default=1)
    params.add_argument("--reference_genome", type=str, help="Reference_genome", choices=["GRCh37", "GRCh38"],
                        default="GRCh38")
    params.add_argument("--test", action='store_true', help="run test with just chrom 22?")
    params.add_argument('--overwrite_checkpoints', type=bool, default=True,
                        help='Overwrite previous pipeline checkpoints?')
    params.add_argument('--run_king', action='store_true', help='Pause pipeline to run King relatedness calculations?')
    params.add_argument('--pc_num', default=10, type=int, help="Number of PCs to calculate.")
    params.add_argument('--verbosity', type=int, default=1,
                        help='Verbosity? Does counts, takes extra time for processing. 0 is least verbose.')
    params.add_argument('--cluster_name', type=str, help='Name of cluster for scaling in pipeline.')
    params.add_argument('--num_preemptible_workers', type=int, default=20,
                        help='Number of preemptible workers for scaling in applicable steps.')

    # Pipeline inputs #
    inputs = parser.add_argument_group("Pipeline inputs and information.")
    inputs.add_argument("-mt", type=str, help="Name of matrix table to run QC pipeline on.")
    inputs.add_argument("--out_name", required=True, type=str, help="Output name ")
    inputs.add_argument("--region", default='europe-west1', help='Region name for scaling in pipeline.')
    inputs.add_argument("--out_dir", type=str, help="Directory to write output data to.")
    inputs.add_argument("--log_dir", type=str, help="Directory to write logs to.")
    inputs.add_argument("--scripts_dir", type=str, help="Directory containing python scripts.")
    inputs.add_argument("--bam_metadata", type=str,
                        help="File containing bam metadata information, if not in sample annotation file.")
    inputs.add_argument("--bam_sample_col", type=str, help="Sample column name in bam metadata file")
    inputs.add_argument("--bam_delim", type=str, default="\t", help="Delimiter in bam metadata file")
    inputs.add_argument("--bam_miss", type=str, default="n/a", help="Missing data string in bam metadata file.")
    inputs.add_argument("--samples_annotation_files", type=str,
                        help="Files to annotate the samples with, comma separated.")
    inputs.add_argument("--samples_col", type=str,
                        help="Name of samples column in sample annotation files. Must be the same in all files.")
    inputs.add_argument("--samples_delim", type=str, default="\t",
                        help="Delimiter in sample annotation files. Must be the same in all files. Default tab.")
    inputs.add_argument("--samples_miss", type=str, default="n/a",
                        help="String for missing values in annotation files, e.g. NA. Must be the same in all files.")
    inputs.add_argument("--fam_id", type=str, help="column name corresponding to sample's family ID. Used in kinship.")
    inputs.add_argument("--pat_id", type=str, help="column name corresponding to sample's paternal ID. Used in kinship.")
    inputs.add_argument("--mat_id", type=str, help="column name corresponding to sample's maternal ID. Used in kinship.")

    # Variant QC thresholds #
    var_thresh = parser.add_argument_group("Variant QC thresholds. If not indicated 'final' or 'low pass' in name, "
                                           "thresholds are used for both low pass variant filtering and final variant"
                                           "filtering.")
    var_thresh.add_argument("--low_pass_p_hwe", default=1e-9, type=float, help="Low pass variant QC HWE cutoff")
    var_thresh.add_argument("--final_p_hwe", default=1e-6, type=float, help="Final variant QC HWE cutoff")
    var_thresh.add_argument("--low_pass_min_call_rate", default=0.8, type=float,
                            help="Low pass variant QC min call rate")
    var_thresh.add_argument("--final_min_call_rate", default=0.9, type=float, help="Final variant QC min call rate")
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


    # LD pruning thresholds #
    ld_thresh = parser.add_argument_group("Thresholds for LD pruning.")
    ld_thresh.add_argument("--r2", default=0.2, type=float, help="r2 correlation cutoff for LD pruning.")
    ld_thresh.add_argument("--bp_window_size", default=500000, type=int,
                           help="Sliding window size for LD pruning in bp.")

    # Kinship thresholds #
    kin_thresh = parser.add_argument_group("Kinship thresholds.")
    kin_thresh.add_argument("--ind_maf", default=0.001, type=float,
                            help="Minor allele frequency cutoff for calculating kinship.")
    kin_thresh.add_argument("--relatives_removal_file", type=str,
                            help="File of related individuals to remove, one sample per line.")

    # Pop outlier options #
    pop_opts = parser.add_argument_group("Options for population outlier removal")
    pop_opts.add_argument("--pop_sd_threshold", default=4, type=int,
                          help="Number of standard deviations from mean of PC1 and PC2 on which we mark samples as "
                               "population outliers.")
    pop_opts.add_argument("--pca_plot_annotations", type=str,
                          help="column annotations with which to annotate PCA plots, comma separated if > 1.")

    # Samples removal options #
    samples_removal = parser.add_argument_group("Arbitrary sample removal options.")
    samples_removal.add_argument("--sample_removal_strings", type=str,
                                 help="string or strings that indicate a sample should be filtered out of the dataset."
                                      "If more than one, comma separated.")
    samples_removal.add_argument("--sample_removal_list", type=str,
                                 help="File containing one sample ID per line that should be removed from the dataset.")

    # Samples QC thresholds #
    samples_thresh = parser.add_argument_group("Samples QC thresholds.")
    samples_thresh.add_argument("--sample_call_rate", type=float,
                                help="Minimum genotype call rate per sample. Default none- samples not filtered on "
                                     "call rate.")
    samples_thresh.add_argument("--chimeras_col", required=True, type=str,
                                help="Column in matrix table or annotation files giving sample chimera percentage")
    samples_thresh.add_argument("--chimeras_max", default=0.05, type=float, help="Max % of chimeras allowed for sample")
    samples_thresh.add_argument("--contamination_col", required=True, type=str,
                                help="Column in matrix table or annotation files giving sample contamination percent.")
    samples_thresh.add_argument("--contamination_max", default=0.05, type=float,
                                help="Max % contamination allowed for sample")
    samples_thresh.add_argument("--batch_col_name", type=str,
                                help="Samples annotation in matrix table or annotation giving batch/cohort for "
                                     "stratified samples QC (TiTv, het/homvar, indel ratios, n singletons).")
    samples_thresh.add_argument("--sampleqc_sd_threshold", default=4, type=int,
                                help="Number of standard deviations from mean sample can deviate on heterozygosity.")

    # Impute sex thresholds #
    sex_thresh = parser.add_argument_group("Impute sex thresholds.")
    sex_thresh.add_argument("--female_threshold", default=0.4, type=float, help="F-stat cutoff for defining female")
    sex_thresh.add_argument("--male_threshold", default=0.8, type=float, help="F-stat cutoff for defining male")

    # Filter by phenotype thresholds #
    pheno_thresh = parser.add_argument_group("Thresholds for filtering by phenotype.")
    pheno_thresh.add_argument("--pheno_col", type=str,
                              help="Samples annotation giving phenotype boolean annotation. Note: samples missing"
                                   "phenotype information are ignored in many cases.")
    pheno_thresh.add_argument("--pheno_call_rate", default=0.95, type=float,
                              help="Min call rate for variant, in cases + controls separately.")

    parsed_args = parser.parse_args(arguments)

    return parsed_args


def check_inputs(parsed_args):
    ################################################################
    # Check that relatives removal file given if run_king is false #
    ################################################################
    if (parsed_args.run_king is False) and (parsed_args.relatives_removal_file is None):
        logging.error("Error! If --run_king is false, then give the file with list of relatives to remove with "
                      "--relatives_removal_file")
        exit(1)

    if (parsed_args.samples_annotation_files is not None) and (parsed_args.samples_col is None):
        logging.error("Error! If --samples_annotation_files given, --samples_col pointing to sample column must also"
                      "be given!")
        exit(1)

    if (parsed_args.num_preemptible_workers is not None) and (parsed_args.cluster_name is None):
        logging.info("Warning! If you want to add preemptible workers during the pipeline, give --cluster_name as well."
                     "No preemptible workers added.")

    ################################
    # Check that input files exist #
    ################################
    files = ['bam_metadata', 'relatives_removal_file', 'sample_removal_list', 'mt', 'scripts_dir',
             'samples_annotation_files']

    for f in files:
        f1 = getattr(parsed_args, f)

        if f == 'samples_annotation_files':
            files = getattr(parsed_args, f)
            if files is not None:
                files = files.strip().split(",")
        else:
            files = [f1]

        if files is not None:
            for file in files:
                if file is not None:
                    if file.endswith("/"):
                        file = file.rstrip("/")

                    if f == "mt":
                        file = file + "/metadata.json.gz"
                    if f == "scripts_dir":
                        file = file + "/exome_qc.py"

                    stat_cmd = ['gsutil', '-q', 'stat', file]
                    status = subprocess.call(stat_cmd)

                    if status != 0:
                        logging.error(f"Error! Input file {file} does not exist!")
                        exit(1)

    ##################################################
    # Check that bucket for output directories exist #
    ##################################################
    outputs = ['out_dir', 'log_dir']

    test_f = open('test.txt', 'w')
    test_f.write("testing bucket exists")
    test_f.close()

    for output in outputs:
        out_dir = getattr(parsed_args, output)
        bucket = "/".join(out_dir.split("/")[0:3])

        test_cmd = ['gsutil', 'cp', 'test.txt', bucket]
        success = subprocess.call(test_cmd)

        if success != 0:
            logging.error(f"Error! Bucket {bucket} for output {out_dir} does not exist!")
            exit(1)
        else:
            rm_cmd = ['gsutil', 'rm', os.path.join(bucket, 'test.txt')]
            subprocess.call(rm_cmd)

    ##############################################################
    # Check that cluster name for adding preemptibles is correct #
    ##############################################################
    if parsed_args.cluster_name is not None:
        cmd = shlex.split(f"gcloud dataproc clusters update {parsed_args.cluster_name} --region {parsed_args.region} "
                          f"--num-secondary-workers 1")

        success = subprocess.call(cmd)

        if success != 0:
            logging.error("Error! Cluster name or region given does not exist!")
            exit(1)
        else:
            cmd = shlex.split(f"gcloud dataproc clusters update {parsed_args.cluster_name} --region {parsed_args.region} "
                              f"--num-secondary-workers 0")

            subprocess.call(cmd)


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
    check_inputs(args)

    # Import python scripts to access helper functions #
    scripts = ["helper_scripts.py", "qc_pipeline_functions.py"]
    for script in scripts:
        hl.spark_context().addPyFile(os.path.join(args.scripts_dir, script))

    import helper_scripts as h
    import qc_pipeline_functions as qc

    ####################
    # Configure logger #
    ####################
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day
    log_file = 'exome-qc_' + timestr + '.txt'

    # Create logger
    root = logging.getLogger()
    root.setLevel(logging.INFO)

    # create formatter
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

    ##################################################
    # Pre-define checkpoint names to load data later #
    ##################################################
    args.samples_annotation = {'prefix': '1', 'suffix': 'after_sample_annotations'}
    args.sample_removal = {'prefix': '2', 'suffix': 'after_sample_removal'}
    args.low_pass_variant_qc = {'prefix': '3', 'suffix': 'after_low_pass_variant_qc'}
    args.maf_ld_prune = {'prefix': '4', 'suffix': 'after_maf_ld_prune'}
    args.find_related_individuals = {'prefix': '5', 'suffix': 'after_finding_relatives'}
    args.find_related_individuals_ld_pruned = {'prefix': '5', 'suffix': 'after_finding_relatives_ld_pruned'}
    args.find_pop_outliers = {'prefix': '6', 'suffix': 'after_finding_pop_outliers'}
    args.impute_sex = {'prefix': '7', 'suffix': 'after_imputing_sex'}
    args.samples_qc = {'prefix': '8', 'suffix': 'after_samples_qc'}
    args.final_variant_qc = {'prefix': '9', 'suffix': 'after_final_variant_qc'}
    args.filter_variants_by_phenotype = {'prefix': '10', 'suffix': 'after_filter_var_by_pheno'}
    args.maf_filter_pcs = {'prefix': '11', 'suffix': 'after_maf_filter_pcs'}
    args.ld_prune_pcs = {'prefix': '12', 'suffix': 'after_ld_prune_pcs'}
    args.final_pc_calculation = {'prefix': '13', 'suffix': 'after_final_pc_calculation'}

    ################
    # Run pipeline #
    ################
    args.cpcounter = 1
    args.output_stem = os.path.join(args.out_dir, args.out_name)
    args.checkpoint_folder = os.path.join(args.out_dir, "checkpoint_mts/")
    args.plot_folder = os.path.join(args.out_dir, "plots")
    args.tmp_counter = 1

    # Load in data according to parameters given
    mt = qc.load_data(args)

    # Annotate samples
    mt = qc.annotate_samples(mt, args)
    # actually exist in the dataset after annotating samples.

    # Samples removal
    mt = qc.remove_samples(mt, args)

    # Low-pass variant QC
    mt = qc.low_pass_var_qc(mt, args)

    # LD Prune and MAF filter dataset for relatedness
    mt, md_ldpruned = qc.maf_ldprune_relatedness(mt, args)

    # Export data to find related individuals in King, if necessary
    mt, mt_ldpruned = qc.find_related_individuals(mt, md_ldpruned, args)

    # Find population outliers (excludes relatives)
    mt = qc.find_pop_outliers(mt, mt_ldpruned, args)

    # Impute sex
    mt = qc.impute_sex(mt, args)

    # Samples QC
    mt = qc.samples_qc(mt, args)

    # Variant QC filtering
    mt = qc.final_variant_qc(mt, args)

    # Filter variants missing by pheno
    mt = qc.find_failing_variants_by_pheno(mt, args)

    # Calculate final PCS
    mt = qc.calculate_final_pcs(mt, args)

    # Export matrix table columns for optmatch
    if not args.run_king:
        mtcols = mt.cols()
        mtcols.export(args.output_stem + '_final_dataset_cols_passingQC.tsv')

    # Send logs and finish-up notice
    logging.info('Pipeline ran successfully! Copying logs and shutting down cluster in 10 minutes.')
    h.copy_logs_output(args.log_dir, log_file=log_file, plot_dir=args.plot_folder)

