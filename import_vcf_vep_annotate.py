"""
Script to import VCF file to Hail format, VEP annotate the dataset, and then save it as a Hail matrix table.

Author: Lea Urpa, August 2020
"""
import os
import argparse


if __name__ == "__main__":
    print('Beginning import pipeline')
    import hail as hl
    import logging

    ##################################
    # Parse arguments for imput data #
    ##################################
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=str, required=True,
                        help="Name of VCF file (or files) to import, comma separated if > 1 file.")
    parser.add_argument("--out_file", type=str, help="Name of matrix table to output.")
    parser.add_argument("--scripts_dir", required=True, type=str,
                        help="Bucket and folder containing the scripts for this library.")
    parser.add_argument("--log_dir", type=str, required=True, help="Bucket and folder where logs should be written to.")
    parser.add_argument("--data_dir", type=str, required=True, help="Bucket and folder where VCFs to import exist.")
    parser.add_argument("--vep_config", default="gs://hail-us-vep/vep85-loftee-gcloud.json",
                        help="Location of Hail VEP configuration json file. Default for cluster started with --vep")
    parser.add_argument("--force_bgz", action='store_true', help="Force blog gzip import? Default true.")
    parser.add_argument("--call_fields", default="PGT", help="Name of genotype call field in VCF, default PGT.")
    parser.add_argument("--test", action='store_true', help="Filters data to just chr 22 for testing purposes.")

    args = parser.parse_args()

    ##################
    # Import scripts #
    ##################
    hl.init()

    scripts = ["helper_scripts.py"]

    for script in scripts:
        hl.spark_context().addPyFile(args.scripts_dir + script)

    import helper_scripts as h

    #####################################
    # Configure logging, define outputs #
    #####################################
    logstem = 'import_vep_annotate-'
    datestr, timestr, log_file = h.configure_logging(logstem=logstem)

    log_dir = os.path.join(args.log_dir, logstem + datestr)

    ####################
    # Import VCF files #
    ####################
    vcf_files = args.vcf.strip().split(",")
    if (len(vcf_files) > 1) and (args.out_file is None):
        logging.error("Error! Must give matrix table file name with --out_file if importing more than one VCF.")
        exit(1)

    logging.info('Importing genotype vcfs.')
    matrix_tables = []

    for vcf in vcf_files:
        vcf_name = os.path.join(args.data_dir, vcf)
        vcf_stem = vcf.replace(".vcf", "")
        vcf_stem = vcf_stem.replace(".gz", "")
        vcf_stem = vcf_stem.replace(".bgz", "")
        mt_name = os.path.join(args.data_dir, vcf_stem)

        hl.import_vcf(vcf_name, force_bgz=args.force_bgz, call_fields=args.call_fields).write(mt_name, overwrite=True)
        mt = hl.read_matrix_table(mt_name)  # faster to write naive mt and then upload again

        if args.test:
            logging.info('Test flag given, filtering to on chrom 22.')
            mt = mt.filter_rows(mt.locus.contig == "22")

        mt = mt.annotate_cols(input_file=vcf)

        logging.info('%s imported count: %s' % (vcf, mt.count()))
        matrix_tables.append(mt)

    ###############################################################
    # If more than 1 VCF given, iteratively combine by union_cols #
    ###############################################################
    if len(matrix_tables) > 1:
        logging.info('Joining matrix tables.')

        for i in range(len(matrix_tables)):
            if i != len(matrix_tables):  # otherwise indexing error for last MT
                mt1 = matrix_tables[i]
                next_mt = matrix_tables[i+1]

                if i == 0:  # if first, combing first two MTs
                    mt = mt1.union_cols(next_mt)
                else:  # Else combine the combined MT with the next MT
                    mt = mt.union_cols(next_mt)

        logging.info('Joined count: ' + str(mt.count()))

    logging.info('Splitting multiallelic variants')
    mt_split = hl.split_multi_hts(mt)
    logging.info('Split count: ' + str(mt_split.count()))

    logging.info('VEP annotating dataset.')
    mt_vep = hl.vep(mt_split, args.vep_location)

    if args.out_file is None:
        out_name = vcf_files[0]
    else:
        out_name = args.out_file

    if args.test:
        out_name = out_name + "_test"

    logging.info('Writing matrix table to bucket.')
    mt_vep.write(os.path.join(args.data_dir, out_name))

    logging.info('Successfully completed import and VEP annotation. Copying logs to bucket and shutting down in 10 min.')
    h.copy_logs_output(log_dir, timestr=timestr, log_file=log_file, out_dir=args.data_dir)