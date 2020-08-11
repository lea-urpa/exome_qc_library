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
    import hail as hl

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
    parser.add_argument("--vep_config", default="gs://hail-us-vep/vep95-GRCh38-loftee-gcloud.json",
                        help="Location of Hail VEP configuration json file. Default for cluster started with --vep")
    parser.add_argument("--reference_genome", default='GRCh37', choices=['GRCh37', 'GRCh38'],
                        help="Reference genome build.")
    parser.add_argument("--chr_prefix", action='store_true', help="Chromosomes are of form 'chr1', NOT '1' etc.")
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

    # Configure logger
    root = logging.getLogger()
    log_formatter = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(filename=log_file, format=log_formatter, level=logging.INFO)

    handler = logging.StreamHandler(sys.stdout)
    root.addHandler(handler)

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
        # Write MT first, then read it from disk #
        vcf_name = os.path.join(args.data_dir, vcf)
        vcf_stem = vcf.replace(".vcf", "")
        vcf_stem = vcf_stem.replace(".gz", "")
        vcf_stem = vcf_stem.replace(".bgz", "") + ".mt"
        mt_name = os.path.join(args.data_dir, vcf_stem)

        # Deal with import if 37 and chroms have chr prefix, or 38 and no prefix.
        if (args.chr_prefix is True) and (args.reference_genome == "GRCh37"):
            recode = {f"chr{i}": f"{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}

        elif (args.chr_prefix is False) and (args.reference_genome == "GRCh38"):
            recode = {f"{i}": f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
        else:
            recode = None

        # If MT does not already exist, load in VCF and then write it to disk
        stat_cmd = ['gsutil', '-q', 'stat', mt_name + "/metadata.json.gz"]
        exists = subprocess.call(stat_cmd)

        if exists == 1:  # stat returns 1 if file/folder does not exist, 0 if it exists
            logging.info(f'Detected mt of input vcf {vcf} does not exist, importing vcf.')
            if recode is None:
                hl.import_vcf(vcf_name, force_bgz=args.force_bgz, call_fields=args.call_fields,
                              reference_genome=args.reference_genome).write(mt_name, overwrite=True)
            else:
                hl.import_vcf(vcf_name, force_bgz=args.force_bgz, call_fields=args.call_fields,
                              reference_genome=args.reference_genome, contig_recoding=recode
                              ).write(mt_name, overwrite=True)
        else:
            logging.info(f"Detected mt of input vcf {vcf} already exists, reading mt directly.")

        mt = hl.read_matrix_table(mt_name)

        if args.test:
            logging.info('Test flag given, filtering to on chrom 22.')
            if args.reference_genome == "GRCh38":
                chrom_code = "chr22"
            else:
                chrom_code = "22"

            mt = mt.filter_rows(mt.locus.contig == chrom_code)

        mt = mt.annotate_cols(input_file=vcf)

        logging.info('%s imported count: %s' % (vcf, mt.count()))
        matrix_tables.append(mt)

    ###############################################################
    # If more than 1 VCF given, iteratively combine by union_cols #
    ###############################################################
    if len(matrix_tables) > 1:
        logging.info('Joining matrix tables.')

        for i in range(len(matrix_tables)):
            if i != len(matrix_tables)-1:  # otherwise indexing error for last MT
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
    mt_vep = hl.vep(mt_split, args.vep_config)

    if args.out_file is None:
        out_name = vcf_files[0]
    else:
        out_name = args.out_file

    if args.test:
        out_name = out_name + "_test"

    logging.info('Writing matrix table to bucket.')
    mt_vep.write(os.path.join(args.data_dir, out_name))

    logging.info('Successfully completed import and VEP annotation. Copying logs to bucket and shutting down in 10 min.')
    h.copy_logs_output(log_dir, timestr=timestr, log_file=log_file, plot_dir=args.data_dir)
