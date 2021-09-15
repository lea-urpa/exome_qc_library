"""
Script to run post-imputation QC, namely info score filtering (optionally overlapping info scores from multiple
files), optionally combining to one VCF file and running batch-wise firth regression.
"""
import logging
import hail as hl

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    files = parser.add_mutually_exclusive_group()
    files.add_argument("--file", default=None, help="Path and file name to the (vcf) dataset.")
    files.add_argument("--file_list", help="Text file containing one input (vcf) file per line", default=None)
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
    args = parser.parse_args()

    hl.init()

    #####################
    # Load data in Hail #
    #####################
    vcf_files = args.vcf.strip().split(",")

    if (len(vcf_files) > 1) and (args.out_file is None):
        logging.error("Error! Must give matrix table file name with --merged_file_name if importing more than one VCF.")
        exit(1)



    #######################################
    # Annotate data sample info, if given #
    #######################################


    ###################################################
    # Find variants not passing info score thresholds #
    ###################################################

    ############################
    # Batch-wise AF comparison #
    ############################