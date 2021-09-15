"""
Script to run post-imputation QC, namely info score filtering (optionally overlapping info scores from multiple
files), optionally combining to one VCF file and running batch-wise firth regression.
"""
import logging
import os
import utils
import hail as hl

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", default=None, help="File name of vcf datasets, comma sep if more than one.")
    parser.add_argument("--data_dir", default="", type=str, help="Directory/folder of VCF file(s).")
    parser.add_argument("--out_dir", type=str, help="Output directory for merged and QCd data.")
    parser.add_argument("--merged_file_name",
                        help="File name of merged dataset if merging more than 1 imputation output.")
    parser.add_argument("--sample_annotation_file", default=None,
                        help="Text file containing ")
    parser.add_argument("--sample_col")
    parser.add_argument("--batch_col")
    parser.add_argument("--info_score_cutoff", default="0.7", type=str,
                        help="Info score (IMPUTE2 info score) threshold, variants below will be removed. "
                             "Can give multiple thresholds, comma separated.")
    parser.add_argument("--reference_genome", default="GRCh38", choices=["GRCh37",  "GRCh38"])
    parser.add_argument("--chr_prefix", action="store_true",
                        help="chromosome codes have 'chr' prefix? Only needed if GRCh37 and chr prefix.")
    parser.add_argument("--force_bgz", default=True, help="Force bgz upload with Hail?")
    parser.add_argument("--call_fields", default="PGT", help="Genotype call field name in VCF")
    parser.add_argument("--test", action="store_true", help="Run test with chromosome 22?")
    parser.add_argument("--force", action="store_true", help="Force re-run through checkpoints?")
    args = parser.parse_args()

    hl.init()

    #####################
    # Load data in Hail #
    #####################
    vcf_files = args.vcf.strip().split(",")

    if (len(vcf_files) > 1) and (args.merged_file_name is None):
        logging.error("Error! Must give matrix table file name with --merged_file_name if importing more than one VCF.")
        exit(1)

    if args.test:
        test_str = "_test"
    else:
        test_str = ""

    if args.merged_file_name is None:
        basename = os.path.basename(vcf_files[0]).replace(".vcf", "").replace(".gz", "").replace(".bgz", "")
    else:
        basename = args.merged_file_name.rstrip("/").replace(".mt", "")

    out_basename = os.path.join(args.out_dir, basename)
    combined_mt_fn = out_basename + f"_combined{test_str}.mt/"

    if (not utils.check_exists(combined_mt_fn)):
        mt = utils.load_vcfs(vcf_files, args.data_dir, force=args.force, test=args.test,
                             chr_prefix=args.chr_prefix, reference_genome=args.reference_genome, force_bgz=args.force_bgz,
                             call_fields=args.call_fields)

        mt = mt.checkpoint(combined_mt_fn, overwrite=True)
    else:
        mt = hl.read_matrix_table(combined_mt_fn)

    logging.info(f"Final matrix table count: {mt.count()}")

    #######################################
    # Annotate data sample info, if given #
    #######################################


    ###################################################
    # Find variants not passing info score thresholds #
    ###################################################

    ############################
    # Batch-wise AF comparison #
    ############################