"""
Script for parsing arguments to exome QC pipeline
Author: Lea Urpa
Date: June 2021
"""
import argparse
import subprocess
import shlex
import logging


def parse_arguments(arguments):
    # TODO think about how we could make this less... extra
    """
    Takes arguments from the command line, runs argparse, and returns args object.

    :param arguments: command line arguments
    :return: args namespace object
    """
    parser = argparse.ArgumentParser(description="Exome sequencing dataset quality control pipeline.")

    # Pipeline parameters #
    params = parser.add_argument_group("Pipeline parameters")
    params.add_argument("--reference_genome", type=str, help="Reference_genome", choices=["GRCh37", "GRCh38"],
                        default="GRCh38")
    params.add_argument("--test", action='store_true', help="run test with just chrom 22?")
    params.add_argument('--run_king', action='store_true', help='Pause pipeline to run King relatedness calculations?')
    params.add_argument('--pc_num', default=10, type=int, help="Number of PCs to calculate.")
    params.add_argument('--verbosity', type=int, default=1,
                        help='Verbosity? Does counts, takes extra time for processing. 0 is least verbose.')
    params.add_argument('--cluster_name', type=str, help='Name of cluster for scaling in pipeline.')
    params.add_argument('--num_secondary_workers', type=int, default=20,
                        help='Number of secondary workers for scaling in applicable steps.')
    params.add_argument('--force', type=bool, default=False, help="Force a re-run of all steps?")

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

    # Add argument for columns to check n samples input files
    columns = ['chimeras_col', 'contamination_col']

    if parsed_args.fam_id is not None:
        columns.append('fam_id')
    if parsed_args.pat_id is not None:
        columns.append('pat_id')
    if parsed_args.mat_id is not None:
        columns.append('mat_id')
    if parsed_args.batch_col_name is not None:
        columns.append('batch_col_name')
    if parsed_args.pheno_col is not None:
        columns.append('pheno_col')

    parsed_args.sample_cols_check = columns

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

        # TODO implement utils check_if_exists function

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

    # TODO write utility function for this
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

    if args.stop_checkpoint is None:
        args.stop_checkpoint = 99