"""
Script to use predicted variant consequence and population information to find putative causal variants for cases.
"""
import time
import os
import logging
import argparse
import sys


def remove_monomorphic(mt, args):
    """
    Takes matrix table and counts the number of non-reference genotypes per variant, removes variants with non-ref GT
    count == 0.
    :param mt: Matrix table to filter
    :param args: arguments for checkpoint output location and name
    :return: returns filtered matrix table
    """
    start0_count = mt.count_rows()
    logging.info(f"Starting number of variants: {start0_count}")
    mt = mt.annotate_rows(non_ref_gt_count=hl.agg.count_where(mt.GT.is_non_ref()))
    mt = mt.filter_rows(mt.non_ref_gt_count > 0, keep=True)

    mt = mt.checkpoint(os.path.join(args.out_dir, args.output_name + "_non_monomorphic_tmp.mt"))

    args.start_count = mt.count_rows()
    logging.info(f"Number of remaining variants after removing monomorphic variants: {args.start_count} "
                 f"({round(args.start_count/start0_count*100, 2)}% of all variants)")

    return mt


def annotate_control_carriers(mt, args):
    """
    Annotate het and hom var carriers for variants in matrix table.

    :param mt: matrix table to annotate/filter
    :param args: arguments for pheno column, etc
    :return: returns annotated matrix table
    """

    logging.info("Annotating het and hom var carrier counts in controls to variants.")

    # Get count of samples that are cases and controls, report
    case_count = mt.aggregate_cols(hl.agg.count_where(mt[args.pheno_col] == True))
    control_count = mt.aggregate_cols(hl.agg.count_where(mt[args.pheno_col] == False))
    missing = mt.aggregate_cols(hl.agg.count_where(~hl.is_defined(mt[args.pheno_col])))

    logging.info(f"Number of controls in dataset: {control_count}")
    logging.info(f"Number of cases in dataset: {case_count}")
    logging.info(f"Samples missing case/control information: {missing}")
    if missing > 0:
        logging.info(f"Warning- samples missing case/control status will be generally ignored in this pipeline.")

    # Annotate control/case het count + homvar count
    mt = mt.annotate_rows(control_het_count=
                          hl.agg.filter((mt[args.pheno_col] == False) & hl.is_defined(mt[args.pheno_col]),
                                        hl.agg.count_where(mt.GT.is_het())))
    mt = mt.annotate_rows(control_homvar_count=
                          hl.agg.filter((mt[args.pheno_col] == False) & hl.is_defined(mt[args.pheno_col]),
                                        hl.agg.count_where(mt.GT.is_hom_var())))
    mt = mt.annotate_rows(case_het_count=
                          hl.agg.filter((mt[args.pheno_col] == True) & hl.is_defined(mt[args.pheno_col]),
                                        hl.agg.count_where(mt.GT.is_het())))
    mt = mt.annotate_rows(case_homvar_count=
                          hl.agg.filter((mt[args.pheno_col] == True) & hl.is_defined(mt[args.pheno_col]),
                                        hl.agg.count_where(mt.GT.is_hom_var())))

    return mt


def annotate_variants(mt, args):
    """
    Annotates matrix table variants with CADD, MPC and gnomad info.
    :param mt: matrix table to annotate
    :param args: arguments for cadd, mpc, and gnomad hail table locations
    :return: returns annotated matrix table
    """
    logging.info("Annotating matrix table with CADD, MPC and Gnomad.")
    # Get variants table from matrix table
    checkpoint_name = os.path.join(args.output_dir, args.output_name + "_annotation_tmp.mt")

    # Annotate variants with CADD
    mt = va.annotate_variants_cadd(mt, args.cadd_ht)

    # Annotate variants with MPC
    mt = va.annotate_variants_mpc(mt, args.mpc_ht)

    # Annotate variants with Gnomad
    mt = va.annotate_variants_gnomad(mt, args.gnomad_ht)

    # Get correct index for gnomad population
    args.gnomad_idx = mt.gnomad_popmax_index_dict.take(1)[0][args.gnomad_population]

    # Annotate variants with Gnomad mismatch
    mt = va.annotate_variants_gnomad_mismatch(mt, args.gnomad_mismatch_ht)

    mt = mt.checkpoint(checkpoint_name, overwrite=True)

    return mt


def putative_causal_geneset(mt, args):
    """
    Looks for LOF and damaging missense variants in a particular geneset, regardless of population frequencies.
    :param mt: matrix table to annotate
    :param args: arguments for output directories, gene set file and name
    :return: returns matrix table, unfiltered but rows annotated for presence in gene set of interest
    """
    logging.info("Searching for putative causal variants in given geneset.")

    # Import gene set if interest, annotate rows with True/False annotation, true if gene name in gene set if interest
    gene_ht = hl.import_table(args.gene_list, no_header=True)
    mt = mt.annotate_rows(**{args.gene_set_name: mt.gene.contains(gene_ht.f0)})

    variants_in_geneset = mt.aggregate_rows(hl.agg.count_where(mt[args.gene_set_name] == True))
    logging.info(f"Number of variants in specified gene set: {variants_in_geneset}")

    # Find LOF and damaging missense variants in these genes
    logging.info(f"Filtering to LOF and damaging missense variants in geneset given in file "
                 f"{os.path.basename(args.gene_list)} with het or hom var cases in this dataset.")
    mt_genes = mt.filter_rows((mt.LOF == True) | ((mt.missense == True) & (mt.MPC > args.mpc_cutoff)) &
                              ((mt.case_het_count > 0) | (mt.case_homvar_count > 0)) &
                              (mt[args.gene_set_name] == True), keep=True)

    geneset_filtered = mt_genes.count_rows()

    logging.info(f"Number of variants in gene set with at least one LOF or damaging missense case carrier: "
                 f"{geneset_filtered} ({round(geneset_filtered/args.start_count*100, 2)}% of non-monomorphic variants)")

    # Get het carriers, by case and control, write matrix table, and pull rows to ht then tsv to write.
    mt_genes = mt_genes.annotate_rows(het_carriers_case=
                                (hl.agg.filter(mt_genes[args.pheno_col] == True) & mt_genes.GT.is_het() &
                                 hl.is_defined(mt_genes.GT), hl.agg.collect(mt_genes.s)))
    mt_genes = mt_genes.annotate_rows(het_carriers_control=
                                  (hl.agg.filter(mt_genes[args.pheno_col] == False) & mt_genes.GT.is_het() &
                                   hl.is_defined(mt_genes.GT), hl.agg.collect(mt_genes.s)))

    # Export results
    mt_genes.write(os.path.join(args.output_dir, args.out_name + "_gene_set_results.mt"))
    gene_variants = mt_genes.rows()

    gene_variants.export(os.path.join(args.output_dir, args.output_name + f"_{args.gene_set_name}_results.tsv"))

    return mt


def putative_causal_dominant(mt, args):
    logging.info("Searching for putative causal variants by dominant inheritance.")

    logging.info(f"Filtering for variants with \n"
                 f"1) if gnomad filters empty, gnomad mismatch between genomes and exomes is False, and gnomad popmax "
                 f"defined, gnomad AC <= {args.max_allowed_carrier_dominant} and no homozygotes in Gnomad\n"
                 f"3) not homvar in any controls in this dataset\n"
                 f"4) het controls in this dataset <= {args.max_allowed_carrier_dominant}\n"
                 f"5) heterozygous in at least one case in this dataset\n"
                 f"6) LOF or missense damaginge (MPC > {args.mpc_cutoff})")

    # Filter out variants with valid gnomad data (filters empty, popmax defined, not mismatch) and AC > threshold and
    # homozygote count > 0
    mt_dom = mt.filter_rows((hl.len(mt.gnomad_filters) == 0) & hl.is_defined(mt.gnomad_popmax[args.gnomad_idx]) &
                            (mt.gnomad_mismatch == False) &
                            (mt.gnomad_popmax.AC[args.gnomad_idx] > args.max_allowed_carrier_dominant) &
                            (mt.gnomad_popmax.homozygote_count[args.gnomad_idx] > 0), keep=False)

    gnomad_filt_count = mt_dom.count_rows()
    logging.info(f"Number of variants after filtering on gnomad annotations: {gnomad_filt_count} "
                 f"({round((args.start_count - gnomad_filt_count)/args.start_count*100, 2)}% "
                 f"of non-monomorphic variants filtered out)")

    # Keep variants with control het count < threshold, controm homvar count 0, case het count > 0, and LOF or
    # damaging missense
    mt_dom = mt_dom.filter_rows((mt_dom.control_het_count <= args.max_allowed_carrier_dominant) &
                                (mt_dom.control_homvar_count == 0) & (mt_dom.case_het_count > 0) &
                                ((mt_dom.LOF == True) | ((mt_dom.missense == True) & (mt_dom.MPC > args.mpc_cutoff))),
                                keep=True)
    putative_causal_dom_count = mt_dom.count_rows()
    logging.info(f"Number of putative causal variants by dominant model in all genes: {putative_causal_dom_count}")
    logging.info(f"{round(putative_causal_dom_count/args.start_count*100, 2)}% of non-monomorphic variants.")

    # Get het carriers, by case and control, write matrix table, and pull rows to ht then tsv to write.
    if putative_causal_dom_count > 0:
        mt_dom = mt_dom.annotate_rows(het_carriers_case=
                                      (hl.agg.filter(mt_dom[args.pheno_col] == True) & mt_dom.GT.is_het() &
                                       hl.is_defined(mt_dom.GT), hl.agg.collect(mt_dom.s)))
        mt_dom = mt_dom.annotate_rows(het_carriers_control=
                                      (hl.agg.filter(mt_dom[args.pheno_col] == False) & mt_dom.GT.is_het() &
                                       hl.is_defined(mt_dom.GT), hl.agg.collect(mt_dom.s)))

        # Export results
        mt_dom.write(os.path.join(args.output_dir, args.out_name + "_dominant_analysis_results.mt"))

        dominant_variants = mt_dom.rows()
        dominant_variants.export(os.path.join(args.output_dir, args.out_name + "_dominant_analysis_results.tsv"))


def putative_causal_recessive(mt, args):
    logging.info("Searching for putative causal variants by recessive inheritance.")

    logging.info(f"Filtering out variants with \n"
                 f"1) if gnomad filters empty, gnomad mismatch between genomes and exomes is False, and gnomad popmax "
                 f"defined, gnomad homozygote count in {args.gnomad_population} < "
                 f"{args.max_allowed_homozygotes_recessive}\n"
                 f"2) and gnomad popmax AF in {args.gnomad_population} < {args.gnomad_AF_cutoff_recessive}\n"
                 f"3) homozygote control count in this dataset < {args.max_allowed_homozygotes_recessive}\n"
                 f"4) LOF or missense damaging (CADD < {args.cadd_cutoff})")

    # Filter out variants with valid gnomad data (filters empty, popmax defined, not mismatch) and
    # homozygote count/AF > threshold
    mt_rec = mt.filter_rows((hl.len(mt.gnomad_filters) == 0) & hl.is_defined(mt.gnomad_popmax[args.gnomad_idx]) &
                            (mt.gnomad_mismatch == False) &
                            (mt.gnomad_popmax.homozygote_count[args.gnomad_idx] > args.max_allowed_homozygotes_recessive)
                            & (mt.gnomad_popmax.AF[args.gnomad_idx] > args.gnomad_AF_cutoff_recessive), keep=False)

    gnomad_filt_count = mt_rec.count_rows()
    logging.info(f"Number of variants after filtering on gnomad annotations: {gnomad_filt_count} "
                 f"({round((args.start_count - gnomad_filt_count)/args.start_count*100, 2)}% of non-monomorphic "
                 f"variants filtered out)")

    # Keep variants with homozygote control count in this dataset < threshold and at least 1 homvar case and LOF or
    # damaging missense
    mt_rec = mt.filter_rows((mt_rec.control_homvar_count <= {args.max_allowed_homozygotes_recessive}) &
                            (mt.case_homvar_count > 0) &
                            ((mt_rec.LOF == True) | ((mt_rec.missense == True) & (mt.CADD_PHRED > args.cadd_cutoff))),
                            keep=True)

    putative_causal_rec_count = mt_rec.count_rows()
    logging.info(f"Number of putative causal variants by recessive model in all genes: {putative_causal_rec_count}")
    logging.info(f"{round(putative_causal_rec_count/args.start_count*100, 2)}% of all non-monomorphic variants.")

    # Get het and homvar carriers, by case and control, write matrix table, and pull rows to ht then tsv to write
    if putative_causal_rec_count > 0:
        mt_rec = mt_rec.annotate_rows(het_carriers_case=
                                      (hl.agg.filter(mt_rec[args.pheno_col] == True) & mt_rec.GT.is_het() &
                                       hl.is_defined(mt_rec.GT), hl.agg.collect(mt_rec.s)))
        mt_rec = mt_rec.annotate_rows(het_carriers_control=
                                      (hl.agg.filter(mt_rec[args.pheno_col] == False) & mt_rec.GT.is_het() &
                                       hl.is_defined(mt_rec.GT), hl.agg.collect(mt_rec.s)))

        mt_rec = mt_rec.annotate_rows(homvar_carriers_case=
                                      (hl.agg.filter(mt_rec[args.pheno_col] == True) & mt_rec.GT.is_hom_var() &
                                       hl.is_defined(mt_rec.GT), hl.agg.collect(mt_rec.s)))
        mt_rec = mt_rec.annotate_rows(homvar_carriers_control=
                                      (hl.agg.filter(mt_rec[args.pheno_col] == False) & mt_rec.GT.is_hom_var() &
                                       hl.is_defined(mt_rec.GT), hl.agg.collect(mt_rec.s)))

        # Export results
        mt_rec.write(os.path.join(args.output_dir, args.out_name + "_recessive_analysis_results.mt"))

        recessive_variants = mt_rec.rows()
        recessive_variants.export(os.path.join(args.output_dir, args.out_name + "_recessive_analysis_results.tsv"))


if __name__ == "__main__":
    import hail as hl
    hl.init()

    ###################
    # Parse arguments #
    ###################
    parser = argparse.ArgumentParser(description="v9 exome sequencing dataset quality control pipeline.")
    parser.add_argument("-mt", type=str, help="Input matrix table to run analysis on")
    parser.add_argument("--pheno_col", required=True, type=str, help="Col annotation giving case status, true/false.")
    parser.add_argument("--mpc_cutoff", type=float, default=2,
                        help="Threshold for damaging missense variants by MPC score, for dominant model.")
    parser.add_argument("--cadd_cutoff", type=float, default=20,
                        help="Threshold for damaging missense variants by CADD score, for recessive model.")
    parser.add_argument("--max_allowed_carrier_dominant", type=int, default=5,
                        help="Maximum allowed number of carriers for a variant for dominant inheritance model.")
    parser.add_argument("--max_allowed_homozygotes_recessive", type=int, default=5, help="Maximum allowed number of "
                        "homozygous carriers for a variant in recessive inheritance model.")
    parser.add_argument("--gnomad_AF_cutoff_recessive", type=float, default=0.1,
                        help="Max allele frequency in Gnomad to consider a variant damaging in recessive model.")
    parser.add_argument("--cadd_ht", required=True, type=str, help="File name of CADD hail table.")
    parser.add_argument("--mpc_ht", required=True, type=str, help="File name of MPC hail table.")
    parser.add_argument("--gnomad_ht", required=True, type=str, help="File name of gnomad hail table.")
    parser.add_argument("--gnomad_population", type=str, default="controls",
                        choices=['gnomad', 'non_topmed', 'non_neuro', 'non_cancer', 'controls'])
    parser.add_argument("--output_name", required=True, type=str, help="Output name for files.")
    parser.add_argument("--output_dir", required=True, type=str, help="Output directory for output files.")
    parser.add_argument("--log_dir", required=True, type=str, help="Output directory for logs.")
    parser.add_argument("--gene_list", type=str,
                        help="Name of file containing one gene name per line, for genes of interest.")
    parser.add_argument("--gene_set_name", default="gene_of_interest")
    args = parser.parse_args()

    scripts = ["variant_annotation.py", "helper_scripts.py"]
    for script in scripts:
        hl.spark_context().addPyFile(os.path.join(args.scripts_dir, script))

    import variant_annotation as va
    import helper_scripts as h

    ####################
    # Configure logger #
    ####################
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day
    args.log_file = 'find-putative-causal-variants_' + timestr + '.txt'

    # Create logger
    root = logging.getLogger()
    root.setLevel(logging.INFO)

    # create formatter
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

    ########################################################
    # Read in matrix table and remove monomorphic variants #
    ########################################################
    full_mt = hl.read_matrix_table(args.mt)
    var_mt = remove_monomorphic(full_mt, args)

    ################################################
    # Annotate rows with control het/homvar counts #
    ################################################
    var_mt = annotate_control_carriers(var_mt, args)

    ################################################
    # Annotate variants table with annotation data #
    ################################################
    var_mt = annotate_variants(var_mt, args)

    ##########################################
    # Run analysis to find putative variants #
    ##########################################
    var_mt = putative_causal_geneset(var_mt, args)
    putative_causal_dominant(var_mt, args)
    putative_causal_recessive(var_mt, args)

    ###########################
    # Copy logs and shut down #
    ###########################
    logging.info('Pipeline ran successfully! Copying logs and shutting down cluster in 10 minutes.')
    h.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.log_dir)
