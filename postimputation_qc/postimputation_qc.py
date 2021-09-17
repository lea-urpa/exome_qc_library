"""
Script to run post-imputation QC, namely info score filtering (optionally overlapping info scores from multiple
files), optionally combining to one VCF file and running batch-wise firth regression.
"""
import logging
import os
import time
import sys
import argparse
from bokeh.io import output_file, save
import samples_annotation as sa
import utils
import hail as hl

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", default=None, help="File name of vcf datasets, comma sep if more than one.")
    parser.add_argument("--data_dir", default="", type=str, help="Directory/folder of VCF file(s).")
    parser.add_argument("--out_dir", type=str, help="Output directory for merged and QCd data.")
    parser.add_argument("--log_dir", type=str, help="Output directory for logs.")
    parser.add_argument("--merged_file_name",
                        help="File name of merged dataset if merging more than 1 imputation output.")
    parser.add_argument("--samples_annotation_files", default=None,
                        help="Text file(s) containing sample information. Comma sep if more than one.")
    parser.add_argument("--info_score_names", type=str, required=True)
    parser.add_argument("--info_score_cutoff", default=0.7, type=float,
                        help="Info score (IMPUTE2 info score) threshold, variants below will be removed. "
                             "Can give multiple thresholds, comma separated.")
    parser.add_argument("--reference_genome", default="GRCh38", choices=["GRCh37",  "GRCh38"])
    parser.add_argument("--chr_prefix", action="store_true",
                        help="chromosome codes have 'chr' prefix? Only needed if GRCh37 and chr prefix.")
    parser.add_argument("--force_bgz", action="store_true", help="Force bgz upload with Hail?")
    parser.add_argument("--call_fields", default="PGT", help="Genotype call field name in VCF")
    parser.add_argument("--test", action="store_true", help="Run test with chromosome 22?")
    parser.add_argument("--force", action="store_true", help="Force re-run through checkpoints?")
    args = parser.parse_args()

    vcf_files = args.vcf.strip().split(",")
    info_score_names = args.info_score_names.strip().split(",")

    if not ((len(info_score_names) == len(vcf_files)) | (len(info_score_names) == 1)) :
        logging.error("--info_score_names must be the same length as the number of VCFs, or length 1")

    ############################
    # Set up logger, init hail #
    ############################
    hl.init()

    ## Configure logger ##
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day
    args.log_file = 'postimputation-qc_' + timestr + '.txt'

    root = logging.getLogger()  # creates logger
    root.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # Add file handler
    fh = logging.FileHandler(args.log_file)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add streaming handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    ## Configure inputs ##
    args.checkpoint_folder = os.path.join(args.out_dir, "checkpoint_mts/")
    args.plot_folder = os.path.join(args.out_dir, "plots")
    stepcount = 1
    if args.test:
        args.test_str = "_test"
    else:
        args.test_str = ""

    #####################
    # Load data in Hail #
    #####################
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

    if (not utils.check_exists(combined_mt_fn)) or args.force:
        mt = utils.load_vcfs(vcf_files, args.data_dir, args.out_dir, force=args.force, test=args.test,
                             chr_prefix=args.chr_prefix, reference_genome=args.reference_genome, force_bgz=args.force_bgz,
                             call_fields=args.call_fields, save_row_annots=True)

        logging.info(f"Writing checkpoint after importing and merging VCFs")
        mt = mt.checkpoint(combined_mt_fn, overwrite=True)
        utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
    else:
        logging.info(f"Detected VCFs imported and uploaded, loading file {combined_mt_fn}")
        mt = hl.read_matrix_table(combined_mt_fn)

    logging.info(f"Final matrix table count: {mt.count()}")

    #######################################
    # Annotate data sample info, if given #
    #######################################
    samples_annotated = out_basename + f"_samples_annotated{test_str}.mt/"

    if args.samples_annotation_files is not None:

        logging.info("Annotating samples with additional information.")
        if (not utils.check_exists(samples_annotated)) or args.force:
            annotation_files = args.samples_annotation_files.strip().split(",")

            for file in annotation_files:
                mt = sa.annotate_cols_from_file(mt, file, args.samples_delim, args.samples_col, args.samples_miss)

            logging.info(f"Writing checkpoint after annotating samples")
            mt = mt.checkpoint(samples_annotated, overwrite=True)
            utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
        else:
            logging.info(f"Detected annotated mt exists, loading file {samples_annotated}")
            mt = mt.read_matrix_table(samples_annotated)

    ###################################################
    # Find variants not passing info score thresholds #
    ###################################################
    logging.info(f"Finding variants with INFO score > {args.info_score_cutoff} in all input chip VCFs.")
    # Given list of info scores, find which structure name they correspond to
    info_score_structs = []
    for name in info_score_names:
        fail = True
        for i in range(len(vcf_files)):
            if i == 0:
                struct_name = "info"
            else:
                struct_name = f"info_{i}"

            try:
                test = mt[struct_name][name]
                info_score_structs.append(f"{struct_name}.{name}")
                fail = False
            except:
                pass
        if fail == True:
            logging.warning(f"Info name {name} not found! This info score field will not be included in calculating "
                            f"whether a variant passes INFO score thresholds in all chip sets.")

    if len(info_score_structs) == 0:
        logging.error("No info score names given were found in the dataset. Check the VCF headers with zcat vcf_name.vcf "
                      "| less for the names, and check your spelling.")
        exit()

    # Pull rows and calculate which variants pass INFO score threshold in all datasets
    row_info = mt.rows()
    row_info = row_info.flatten()

    row_info = row_info.annotate(
        passing_all_info=hl.all([(row_info[x][0] > args.info_score_cutoff) for x in info_score_structs])
    )

    # Plot hists for info scores for each chip set
    for name in info_score_structs:
        output_file(f"{datestr}_info_score_hist_{name}")
        info_hist = row_info.aggregate(hl.expr.aggregators.hist(row_info[name][0], 0, 1, 50))
        p = hl.plot.histogram(info_hist, legend='IMPUTE2 Info Score', title=f'INFO scores chipset {name}')
        save(p)

    info_count = row_info.aggregate(hl.agg.counter(row_info.passing_all_info))
    logging.info(f"Number of variants passing/failing INFO thresholds at all chipsets:\n"
                 f"passing:{info_count[True]}, failing:{info_count[False]}")

    row_info = row_info.key_by("locus", "alleles")
    mt = mt.annotate_rows(passing_all_info=row_info[mt.locus, mt.alleles].passing_all_info)


    ############################
    # Batch-wise AF comparison #
    ############################
    logging.info("Running pairwise AF comparison (Fisher or Chisq) for each variant, pairwise between input chips.")
    # Annotate ref and alt allele count for all input chip sets
    input_files_short = []
    for input_file in vcf_files:
        annot_name = input_file.replace("/", "").replace("*", "").replace(".vcf", "").replace(".gz", "")
        input_files_short.append(annot_name)
        mt = mt.annotate_rows(
            **{f"{annot_name}_AC_ref": hl.agg.filter(
                mt.input_file == input_file,
                hl.int32((hl.agg.count_where(mt.GT.is_hom_ref()) * 2) + hl.agg.count_where(mt.GT.is_het())))}
        )
        mt = mt.annotate_rows(
            **{f"{annot_name}_AC_alt": hl.agg.filter(
                mt.input_file == input_file,
                hl.int32((hl.agg.count_where(mt.GT.is_hom_var()) * 2) + hl.agg.count_where(mt.GT.is_het())))}
        )

    # Run pairwise AF test for each variant
    for pair in utils.get_upper_triangle(input_files_short):
        chip1 = pair[0]
        chip2 = pair[1]
        mt = mt.annotate_rows(
            **{f"af_test_{chip1}_{chip2}":
                   hl.contingency_table_test(mt[f"{chip1}_AC_ref"], mt[f"{chip2}_AC_ref"],
                                             mt[f"{chip1}_AC_alt"], mt[f"{chip2}_AC_alt"],
                                             min_cell_count=500)
               }
        )

        output_file(f"{datestr}_AF_comparison_pvalue_hist_{chip1}_{chip2}")
        info_hist = row_info.aggregate(hl.expr.aggregators.hist(mt[f"af_test_{chip1}_{chip2}"], 0, 1, 50))
        p = hl.plot.histogram(info_hist, legend='p value',
                              title=f'Distribution of p values for AF difference test \nbetween {chip1} and {chip2}')
        save(p)

    final_checkpoint = out_basename + f"_info_af_comparison_final{test_str}.mt/"

    mt = mt.checkpoint(final_checkpoint)
    utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)
    logging.info("Pipeline completed successfully!")