"""
Script to import VCF file to Hail format, VEP annotate the dataset, and then save it as a Hail matrix table.

Author: Lea Urpa, August 2020
"""
# TODO mark all temp files with tmp in output name. Mark final output file

if __name__ == "__main__":
    import os
    import argparse
    import sys
    import subprocess
    import logging
    import time
    import hail as hl
    import utils

    ##################################
    # Parse arguments for imput data #
    ##################################
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vcf", type=str, required=True,
        help="File name of vcf datasets, comma sep if more than one. NOTE: wildcard input only works with VCF files "
             "with identical sample IDs, split by variants, unless you add --sample_split flag.")
    parser.add_argument(
        "--sample_split", action="store_true",
        help="Indicates wildcard in --vcf shows sample files, with each VCF containing one sample ID.")
    parser.add_argument("--out_file", type=str, help="Name of matrix table to output.")
    parser.add_argument("--region", default="europe-west1", help="Name of region for checking dataproc in correct region.")
    parser.add_argument("--project", required=True, help="Project name for requester pays config.")
    parser.add_argument("--log_dir", type=str, required=True, help="Location where logs should be written to.")
    parser.add_argument("--log_debug", action="store_true", help="Print debugging information?")
    parser.add_argument("--data_dir", type=str, required=True, help="Location where VCFs to import exist.")
    parser.add_argument("--out_dir", type=str, required=True, help="Location to write combined + vep annotated mt")
    parser.add_argument("--reference_genome", default='GRCh37', choices=['GRCh37', 'GRCh38'],
                        help="Reference genome build.")
    parser.add_argument("--chr_prefix", action='store_true', help="Chromosomes are of form 'chr1', NOT '1' etc.")
    parser.add_argument("--force_bgz", action='store_true', help="Force blog gzip import?")
    parser.add_argument("--force_load", action='store_true', help="Force loading of gzipped vcf?")
    parser.add_argument("--call_fields", default="PGT", help="Name of genotype call field in VCF, default PGT.")
    parser.add_argument("--test", action='store_true', help="Filters data to just chr 22 for testing purposes.")
    parser.add_argument("--force", action='store_true', help="Force re-run of all steps?")
    parser.add_argument("--liftover_37_to_38", action='store_true', help="Liftover GRCh37 data to GRCh38 before VEP?")

    args = parser.parse_args()

    if 'europe' in args.region:
        vep_bucket = "hail-eu-vep"
        args.vep_config = "gs://hail-eu-vep/vep95-GRCh38-loftee-gcloud.json"
    elif 'us' in args.region:
        vep_bucket = "hail-us-vep"
        args.vep_config = "gs://hail-us-vep/vep95-GRCh38-loftee-gcloud.json"
    else:
        logging.error("Error- region not in Europe or US. Are you sure you want to VEP annotate with files from "
                      "another region? Network egress charges add up VERY quickly. Exiting.")
        exit()

    hl.init(spark_conf={
        'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
        'spark.hadoop.fs.gs.requester.pays.buckets': vep_bucket,
        'spark.hadoop.fs.gs.requester.pays.project.id': args.project
        })

    #####################################
    # Configure logging, define outputs #
    #####################################
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day
    log_file = 'import_vep_annotate-' + timestr + '.txt'

    root = logging.getLogger()  # creates logger
    root.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Add file handler
    fh = logging.FileHandler(log_file)
    if args.log_debug:
        fh.setLevel(logging.DEBUG)
    else:
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
    for file_url in [args.log_dir, args.out_dir, args.data_dir]:
        utils.check_regions(args.region, file_url)

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
        mt = utils.load_vcfs(vcf_files, args.data_dir, args.out_dir, combined_mt_fn, force=args.force, test=args.test,
                             chr_prefix=args.chr_prefix, reference_genome=args.reference_genome, force_bgz=args.force_bgz,
                             call_fields=args.call_fields, force_load=args.force_load, sample_split=args.sample_split)

    else:
        mt = hl.read_matrix_table(combined_mt_fn)

    logging.info(f"Final matrix table count: {mt.count()}")

    ###############################
    # Split multiallelic variants #
    ###############################
    split_fn = out_basename + f"_split{test_str}.mt/"

    if (not utils.check_exists(split_fn)) or args.force:
        logging.info('Splitting multiallelic variants')
        mt_split = hl.split_multi_hts(mt)
        mt_split = mt_split.checkpoint(split_fn, overwrite=True)
    else:
        logging.info("Detected split mt exists, loading that.")
        mt_split = hl.read_matrix_table(split_fn)
    logging.info('Split count: ' + str(mt_split.count()))

    ################################################
    # If indicated, liftover from GRCh37 to GRCh38 #
    ################################################
    if args.liftover_37_to_38:
        lifted_fn = f"{out_basename}_lifted_to_GRCh38{test_str}.mt/"

        if (not utils.check_exists(lifted_fn)) or args.force:
            logging.info(f"Lifting over file from GRCh37 to GRCh38")
            rg37 = hl.get_reference('GRCh37')
            rg38 = hl.get_reference('GRCh38')
            rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

            mt_split = mt_split.annotate_rows(new_locus=hl.liftover(mt_split.locus, 'GRCh38', include_strand=True),
                                              old_locus=mt_split.locus)

            logging.info("Variants that exist in GRCh38:")
            logging.info(mt_split.aggregate_rows(hl.agg.counter(hl.is_defined(mt_split.new_locus))))
            logging.info("Variants that exist in GRCh38 and are not flipped strand:")
            logging.info(mt_split.aggregate_rows(hl.agg.counter(
                hl.is_defined(mt_split.new_locus) & ~mt_split.new_locus.is_negative_strand
            )))
            logging.info("Variants that are strand flipped and do not exist in GRCh38 will be removed.")

            mt_split = mt_split.filter_rows(hl.is_defined(mt_split.new_locus) & ~mt_split.new_locus.is_negative_strand)
            mt_split = mt_split.key_rows_by(locus=mt_split.new_locus.result, alleles=mt_split.alleles)

            mt_split = mt_split.checkpoint(lifted_fn, overwrite=True)

        else:
            logging.info("Detected lifted file already exists, loading that.")
            mt_split = hl.read_matrix_table(lifted_fn)

        args.reference_genome = "GRCh38"
        args.vep_config = "gs://hail-us-vep/vep95-GRCh38-loftee-gcloud.json"

    ########################
    # VEP annotate dataset #
    ########################
    vep_fn = out_basename + f"_vep_annot{test_str}.mt/"

    logging.info('VEP annotating dataset.')
    mt_vep = hl.vep(mt_split, args.vep_config)

    logging.info('Writing matrix table to bucket.')
    mt_vep.write(vep_fn, overwrite=True)

    logging.info('Successfully completed import and VEP annotation. Copying logs to bucket and shutting down in 10 min.')
    utils.copy_logs_output(log_dir, log_file=log_file, plot_dir=args.data_dir)
