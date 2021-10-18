"""
Script to import VCF file to Hail format, VEP annotate the dataset, and then save it as a Hail matrix table.

Author: Lea Urpa, August 2020
"""


if __name__ == "__main__":
    print('Beginning import pipeline')
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
    parser.add_argument("--vcf", type=str, required=True,
                        help="Name of VCF file (or files) to import, comma separated if > 1 file.")
    parser.add_argument("--out_file", type=str, help="Name of matrix table to output.")
    parser.add_argument("--log_dir", type=str, required=True, help="Location where logs should be written to.")
    parser.add_argument("--data_dir", type=str, required=True, help="Location where VCFs to import exist.")
    parser.add_argument("--out_dir", type=str, required=True, help="Location to write combined + vep annotated mt")
    parser.add_argument("--vep_config", default="gs://hail-us-vep/vep95-GRCh38-loftee-gcloud.json",
                        help="Location of Hail VEP configuration json file. Default for cluster started with --vep")
    parser.add_argument("--reference_genome", default='GRCh37', choices=['GRCh37', 'GRCh38'],
                        help="Reference genome build.")
    parser.add_argument("--chr_prefix", action='store_true', help="Chromosomes are of form 'chr1', NOT '1' etc.")
    parser.add_argument("--force_bgz", action='store_true', help="Force blog gzip import? Default true.")
    parser.add_argument("--call_fields", default="PGT", help="Name of genotype call field in VCF, default PGT.")
    parser.add_argument("--test", action='store_true', help="Filters data to just chr 22 for testing purposes.")
    parser.add_argument("--force", action='store_true', help="Force re-run of all steps?")

    args = parser.parse_args()

    hl.init()

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

    if (not utils.check_exists(combined_mt_fn)) or args.force:
        mt = utils.load_vcfs(vcf_files, args.data_dir, args.out_dir, force=args.force, test=args.test,
                             chr_prefix=args.chr_prefix, reference_genome=args.reference_genome, force_bgz=args.force_bgz,
                             call_fields=args.call_fields)

        mt = mt.checkpoint(combined_mt_fn, overwrite=True)
    else:
        mt = hl.read_matrix_table(combined_mt_fn)

    logging.info(f"Final matrix table count: {mt.count()}")

    ##############################################
    # Split multiallelic variants + VEP annotate #
    ##############################################
    split_fn = out_basename + f"_split{test_str}.mt/"
    vep_fn = out_basename + f"_vep_annot{test_str}.mt/"

    if (not utils.check_exists(split_fn)) or args.force:
        logging.info('Splitting multiallelic variants')
        mt_split = hl.split_multi_hts(mt)
        mt_split = mt_split.checkpoint(split_fn, overwrite=True)
    else:
        logging.info("Detected split mt exists, loading that.")
        mt_split = hl.read_matrix_table(split_fn)
    logging.info('Split count: ' + str(mt_split.count()))

    logging.info('VEP annotating dataset.')
    mt_vep = hl.vep(mt_split, args.vep_config)

    logging.info('Writing matrix table to bucket.')
    mt_vep.write(vep_fn, overwrite=True)

    logging.info('Successfully completed import and VEP annotation. Copying logs to bucket and shutting down in 10 min.')
    utils.copy_logs_output(log_dir, log_file=log_file, plot_dir=args.data_dir)
