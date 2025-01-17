"""
Script for doing exome sequencing data quality control from a Hail matrix table input that has been VEP annotated.

Author: Lea Urpa, August 2020
This pipeline is dedicated to Thao and the Get Down Stay Down, whose album A Man Alive was the mojo for writing
"""
import sys
import time
import logging
import os
import hail as hl
from bokeh.io import output_file, save
from parse_arguments import parse_arguments, check_inputs
import utils
import samples_annotation as sa
import variant_qc as vq
import samples_qc as sq
import variant_annotation as va


def initialize_logger(args):
    """Initializes and configures the logger."""
    datestr = time.strftime("%Y.%m.%d")
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")
    args.log_file = f'exome_qc_{timestr}.txt'

    root = logging.getLogger()
    root.setLevel(logging.INFO if not args.log_debug else logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    file_handler = logging.FileHandler(args.log_file)
    file_handler.setFormatter(formatter)
    root.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    root.addHandler(console_handler)


def verify_inputs(args):
    """Verifies that all input arguments are correct, and that the input file GC bucket matches the dataproc cluster"""
    for folder in [args.log_dir, args.out_dir, args.mt]:
        utils.check_regions(args.region, folder)

    check_inputs(args)


def load_and_annotate_samples(args):
    """Loads and annotates samples based on input arguments."""
    checkpoint_path = os.path.join(args.checkpoint_dir, f"samples_annotated.mt")

    if not utils.check_exists(checkpoint_path) or args.force:
        logging.info(f"Loading matrix table: {args.mt}")
        mt = hl.read_matrix_table(args.mt)
        mt.annotate_globals(original_mt_input={'file': args.mt})
        utils.check_vep(mt)
        utils.check_entry(mt)

        if args.test:
            utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)
            mt = utils.create_test_dataset(mt, args.reference_genome, args.mt, args.out_dir)
            utils.remove_secondary(args.cluster_name, args.region)

        if args.samples_annotation_files:
            logging.info('Annotating samples.')
            for annotation_file in args.samples_annotation_files.split(","):
                mt = sa.annotate_cols_from_file(
                    mt,
                    annotation_file,
                    args.samples_delim,
                    args.samples_col,
                    args.samples_miss
                )

        if args.bam_metadata:
            mt = sa.annotate_cols_from_file(
                mt,
                args.bam_metadata,
                args.bam_delim,
                args.bam_sample_col,
                args.bam_miss
            )

        mt = mt.checkpoint(checkpoint_path, overwrite=True)
    else:
        logging.info("Detected sample-annotated mt exists, skipping samples annotation.")
        mt = hl.read_matrix_table(checkpoint_path)

    return mt


def verify_annotations(args):
    """Confirms necessary annotations exist in the data after annotation."""
    for column_name in args.sample_cols_check:
        col = getattr(args, column_name)
        try:
            test = hl.is_defined(mt[col])
        except Exception as e:
            logging.error(f"Error! Given column annotation {col} does not actually exist after inputting sample "
                          f"annotations.")
            logging.error(e)
            exit(1)


def case_specific_variant_qc(mt, args, qc_type):
    """Runs case-specific variant QC."""
    if args.pheno_col is not None:
        logging.info("Running case-status specific variant QC")
        mt = vq.find_variants_failing_by_pheno(mt, ab_allowed_dev_het=args.ab_allowed_dev_het,
                                               pheno_call_rate=args.pheno_call_rate, prefix=qc_type)
    else:
        logging.info("Phenotype column not given, skipping filtering variants by phenotype.")

    return mt


def run_variant_qc(mt, args, qc_type):
    """Performs variant quality control."""
    qc_path = os.path.join(args.out_dir, f"{qc_type}_variant_qc.mt")

    if not utils.check_exists(qc_path) or args.force:
        logging.info(f"Running {qc_type} variant and genotype QC.")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)

        if qc_type == "final":
            sex_aware_call_rate = True
            samples_qc = True
        else:
            sex_aware_call_rate = False
            samples_qc = False

        mt = vq.variant_quality_control(
            mt,
            qc_path,
            annotation_prefix=qc_type,
            min_dp=args.min_dp,
            min_gq=args.min_gq,
            max_het_ref_reads=args.max_het_ref_reads,
            min_het_ref_reads=args.min_het_ref_reads,
            min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
            max_hom_alt_ref_reads=args.max_hom_alt_ref_reads,
            call_rate=args[f"{qc_type}_min_call_rate"],
            p_hwe=args[f"{qc_type}_p_hwe"],
            snp_qd=args.snp_qd,
            indel_qd=args.indel_qd,
            ab_allowed_dev_het=args.ab_allowed_dev_het,
            count_failing=args.count_failing,
            sex_aware_call_rate=sex_aware_call_rate,
            pheno_col=args.pheno_col,
            samples_qc=samples_qc,
            force=args.force
        )

        if qc_type == "final":
            mt = case_specific_variant_qc(mt, args, qc_type)

        mt = mt.checkpoint(qc_path, overwrite=True)

    else:
        logging.info(f"Detected {qc_type} variant and genotype QC already performed.")
        mt = hl.read_matrix_table(qc_path)

    return mt


def plot_variant_stats(mt, qc_type):
    output_file(f"{qc_type}_mean_het_ab_hist.html")
    ab_hist = mt.aggregate_rows(hl.agg.hist(mt[f"{qc_type}_het_ab_stats"].mean, 0, 1, 50))
    p = hl.plot.histogram(ab_hist, legend='het ref read ratio', title='Mean het read ratio per var (passing GTs)')
    save(p)

    output_file(f"{qc_type}_initial_call_rate.html")
    cr_hist_1 = mt.aggregate_rows(hl.agg.hist(mt.low_pass_initial_call_rate, 0, 1, 50))
    p1 = hl.plot.histogram(cr_hist_1, legend='call rate', title="variant call rate, before GT filters")
    save(p1)


def filter_failing_variants_and_genotypes(mt, args, checkpoint_name, qc_type, unfilter_entries, keep_hwe,
                                          filter_samples=False, pop_outliers=False):
    """Filters failing variants and genotypes, and filters out low MAF variants."""
    if not utils.check_exists(checkpoint_name) or args.force:
        logging.info(f"Filtering out failing variants and genotypes.")
        mt = sq.filter_failing(
            mt,
            checkpoint_name,
            prefix=qc_type,
            entries=True,
            variants=True,
            samples=filter_samples,
            unfilter_entries=unfilter_entries,
            pheno_qc=False,
            keep_hwe=keep_hwe,
            min_dp=args.min_dp,
            min_gq=args.min_gq,
            max_het_ref_reads=args.max_het_ref_reads,
            min_het_ref_reads=args.min_het_ref_reads,
            min_hom_ref_ref_reads=args.min_hom_ref_ref_reads,
            max_hom_alt_ref_reads=args.max_hom_alt_ref_reads,
            force=args.force,
            pop_outliers=pop_outliers
        )

        mt = mt.checkpoint(checkpoint_name, overwrite=True)
    else:
        mt = hl.read_matrix_table(checkpoint_name)

    return mt


def filter_maf(mt, args, checkpoint_name, qc_type):
    """Filters variants with MAF below cutoff"""
    if not utils.check_exists(checkpoint_name) or args.force:
        logging.info(f"Filtering to variants with MAF < {args.ind_maf}")
        mt = vq.maf_filter(mt, maf=args.ind_maf, varqc_annot_name=f"{qc_type}_variant_qc")

        mt = mt.checkpoint(checkpoint_name, overwrite=True)
    else:
        mt = hl.read_matrix_table(checkpoint_name)

    return mt


def downsample(mt, args, checkpoint_name):
    """Downsamples MT, first with ld pruning and then randomly downsampling if above target variant count."""
    if not utils.check_exists(checkpoint_name) or args.force:
        mt = vq.downsample_variants(
            mt,
            target_count=80000,
            checkpoint_name=checkpoint_name,
            r2=args.r2,
            bp_window_size=args.bp_window_size,
            ld_prune=True
        )
    else:
        mt = hl.read_matrix_table(checkpoint_name)

    return mt


def filter_to_autosomes(mt, args):
    """Filters dataset to just autosomes"""
    if args.reference_genome == "GRCh38":
        autosomes = ["chr" + str(i) for i in range(1, 23)]
    else:
        autosomes = [str(i) for i in range(1, 23)]

    mt = mt.filter_rows(hl.literal(autosomes).contains(mt.locus.contig))

    return mt


def calculate_relatedness(mt, args, qc_type):
    """Calculates relatedness and annotates the matrix table."""
    relatedness_path = os.path.join(args.out_dir, "relatedness_annotated.mt")
    relatedness_downsampled_path = os.path.join(args.out_dir, "relatedness_downsampled.mt")

    if not utils.check_exists(relatedness_path) or args.force:
        logging.info("Calculating relatedness")
        filtered_failing_path = os.path.join(args.out_dir, "filtered_variants.mt")
        maf_filtered_path = os.path.join(args.out_dir, "maf_filtered.mt")
        downsampled_path = os.path.join(args.out_dir, "downsampled.mt")

        mt_downsampled = filter_failing_variants_and_genotypes(
            mt,
            args,
            checkpoint_name=filtered_failing_path,
            qc_type=qc_type,
            unfilter_entries= True,
            keep_hwe=False,
            filter_samples=False,
            pop_outliers=False
        )
        mt_downsampled = filter_maf(mt_downsampled, args, maf_filtered_path, qc_type)
        mt_downsampled = downsample(mt_downsampled, args, downsampled_path)
        mt_downsampled = filter_to_autosomes(mt_downsampled, args)

        related_to_remove, related_info_ht = sq.king_relatedness(
            mt_downsampled, relatedness_path, kinship_threshold=args.kinship_threshold,
            pheno_col=args.pheno_col, force=args.force,
            cluster_name=args.cluster_name, num_secondary_workers=args.num_secondary_workers,
            region=args.region
        )

        mt = mt.annotate_cols(
            related_to_remove=hl.if_else(hl.literal(related_to_remove).contains(mt.s), True, False),
            related_graph_id=related_info_ht[mt.s].related_graph_id,
            related_num_connections=hl.or_else(related_info_ht[mt.s].related_num_connections, 0)
        )

        mt_downsampled = mt_downsampled.annotate_cols(
            related_to_remove=hl.if_else(hl.literal(related_to_remove).contains(mt_downsampled.s), True, False),
            related_graph_id=related_info_ht[mt_downsampled.s].related_graph_id,
            related_num_connections=hl.or_else(related_info_ht[mt_downsampled.s].related_num_connections, 0)
        )

        mt = mt.checkpoint(relatedness_path, overwrite=True)
        mt_downsampled = mt_downsampled.checkpoint(relatedness_downsampled_path, overwrite=True)
    else:
        logging.info("Detected relatedness already calculated.")
        mt = hl.read_matrix_table(relatedness_path)
        mt_downsampled = hl.read_matrix_table(relatedness_downsampled_path)

    return mt, mt_downsampled


def find_population_outliers(mt, mt_downsampled, args):
    """Identifieds population outliers using the PCA method."""
    pop_outliers_path = os.path.join(args.out_dir, "pop_outliers_annotated.mt")
    pop_outliers_downsampled_path = os.path.join(args.out_dir, "pop_outliers_downsampled.mt")

    if not utils.check_exists(pop_outliers_downsampled_path) or args.force:
        logging.info("Finding population outliers with PCA method.")
        utils.add_secondary(args.cluster_name, args.num_secondary_workers, args.region)
        pop_outliers = sq.find_pop_outliers(
            mt_downsampled,
            pop_outliers_downsampled_path,
            pop_sd_threshold=args.pop_sd_threshold,
            plots=args.pca_plots,
            max_iter=args.max_iter,
            reference_genome=args.reference_genome,
            pca_plot_annotations=args.pca_plot_annotations
        )

        mt = mt.annotate_cols(pop_outlier_sample=hl.if_else(hl.literal(pop_outliers).contains(mt.s), True, False))
        mt_downsampled = mt_downsampled.annotate_cols(
            pop_outlier_sample=hl.if_else(hl.literal(pop_outliers).contains(mt_downsampled.s), True, False))

        mt = mt.checkpoint(pop_outliers_path, overwrite=True)
        mt_downsampled = mt_downsampled.checkpoint(pop_outliers_downsampled_path, overwrite=True)
    else:
        logging.info("Detected population outliers already found.")
        mt = hl.read_matrix_table(pop_outliers_path)
        mt_downsampled = hl.read_matrix_table(pop_outliers_downsampled_path)

    return mt, mt_downsampled


def annotate_variants(mt, args):
    """Adds custom variant annotations from VEP."""
    annotated_path = os.path.join(args.out_dir, "variants_annotated.mt")

    if not utils.check_exists(annotated_path) or args.force:
        logging.info("Adding custom annotations to variants.")
        mt = va.annotate_variants(mt)
        mt = mt.checkpoint(annotated_path, overwrite=True)
    else:
        mt = hl.read_matrix_table(annotated_path)

    return mt


def impute_sex(mt, args, qc_type):
    """Imputes sex of variants."""
    sex_imputed_path = os.path.join(args.out_dir, "sex_imputed.mt")
    filtered_failing_keephwe_path = os.path.join(args.out_dir, "filtered_variants_keephwe_failing.mt")

    if not utils.check_exists(sex_imputed_path) or args.force:
        logging.info("Imputing sex and calculating sex-aware sample annotations.")
        mt_filtered = filter_failing_variants_and_genotypes(
            mt,
            args,
            checkpoint_name=filtered_failing_keephwe_path,
            qc_type=qc_type,
            unfilter_entries=False,
            keep_hwe=True,
            filter_samples=False,
            pop_outliers=False
        )

        imputed_sex = sq.impute_sex_plot(
            mt_filtered,
            female_threshold=args.female_threshold,
            male_threshold=args.male_threshold
        )

        mt_filtered = mt_filtered.annotate_cols(is_female_imputed=imputed_sex[mt_filtered.s].is_female)
        mt = mt.annotate_cols(is_female_imputed=imputed_sex[mt.s].is_female, f_stat=imputed_sex[mt.s].f_stat)

        mt_filtered = sa.sex_aware_sample_annotations(mt_filtered)
        mt = mt.annotate_cols(sexaware_sample_call_rate=mt_filtered.cols()[mt.s].sexaware_sample_call_rate)

        mt = mt.checkpoint(sex_imputed_path, overwrite=True)
    else:
        mt = hl.read_matrix_table(sex_imputed_path)

    return mt


def run_samples_qc(mt, args):
    """Runs quality control on samples."""
    qc_path = os.path.join(args.out_dir, "samples_qc.mt")
    filtered_failing_path = os.path.join(args.out_dir, "filtered_variants.mt")

    if not utils.check_exists(qc_path) or args.force:
        logging.info("Running samples QC.")

        mt_filtered = filter_failing_variants_and_genotypes(
            mt,
            args,
            checkpoint_name=filtered_failing_path,
            qc_type="low_pass",
            unfilter_entries=False,
            keep_hwe=True,
            filter_samples=False,
            pop_outliers=False
        )

        utils.remove_secondary(args.cluster_name, args.region)

        mt = sq.samples_qc(
            mt_filtered,
            mt,
            qc_path,
            count_failing=args.count_failing,
            sample_call_rate=args.sample_call_rate,
            chimeras_col=args.chimeras_col,
            chimeras_max=args.chimeras_max,
            contamination_col=args.contamination_col,
            contamination_max=args.contamination_max,
            batch_col_name=args.batch_col_name,
            sampleqc_sd_threshold=args.sampleqc_sd_threshold,
            pheno_col=args.pheno_col
        )

        mt = mt.checkpoint(qc_path, overwrite=False)
    else:
        logging.info("Detected samples QC run.")
        mt = hl.read_matrix_table(qc_path)

    return mt


def plot_pcs(mt, args):
    """Plot principal components, with or without annotations"""
    if args.pca_plot_annotations is not None:
        try:
            pca_annotations = args.pca_plot_annotations.strip().split(",")
            label_dict = {i: mt[i] for i in pca_annotations}

            output_file(f"final_pcs_plot.html")
            p = hl.plot.scatter(mt.pc1, mt.pc2, label=label_dict, title="Final PCs", collect_all=True)
            save(p)
        except Exception as e:
            logging.error(f"Error! Creating PCA plots with labels failed. Are the label categories you provided"
                          f" really in the data? labels provided: {args.pca_plot_annotations}. Plotting without "
                          f"labels")
            logging.error(e)
            output_file(f"final_pcs_plot.html")
            p = hl.plot.scatter(mt.pc1, mt.pc2, title="Final principal components", collect_all=True)
            save(p)
    else:
        output_file(f"final_pcs_plot.html")
        p = hl.plot.scatter(mt.pc1, mt.pc2, title="Final principal components", collect_all=True)
        save(p)


def calculate_final_pcs(mt, args, qc_type):
    """Calculates final principal components and projects them onto related samples."""
    pcs_annotated_path = os.path.join(args.out_dir, "final_pcs.mt")

    if not utils.check_exists(pcs_path) or args.force:
        logging.info("Calculating principal components for passing samples with passing variants and genotypes.")
        mt_filtered_path = os.path.join(args.out_dir, "final_filtered.mt")

        filtered_failing_path = os.path.join(args.out_dir, "final_passing_variants_samples_only.mt")
        maf_filtered_path = os.path.join(args.out_dir, "maf_filtered_pcs.mt")
        downsampled_path = os.path.join(args.out_dir, "downsampled_pcs.mt")
        pcs_path = os.path.join(args.out_dir, "pcs.mt/")

        mt_downsampled = filter_failing_variants_and_genotypes(
            mt,
            args,
            checkpoint_name=filtered_failing_path,
            qc_type=qc_type,
            unfilter_entries=True,
            keep_hwe=False,
            filter_samples=True,
            pop_outliers=True
        )
        mt_downsampled = filter_maf(mt_downsampled, args, maf_filtered_path, qc_type)
        mt_downsampled = downsample(mt_downsampled, args, downsampled_path)
        mt_downsampled = filter_to_autosomes(mt_downsampled, args)

        scores, related_scores = sq.project_pcs_relateds(
            mt_downsampled, pcs_path, covar_pc_num=args.pc_num, reference_genome=args.reference_genome
        )

        mt = mt.annotate_cols(**{'pc' + str(k + 1): scores[mt.s].scores[k]
                                 for k in range(args.pc_num)})
        mt = mt.annotate_cols(**{'pc' + str(k + 1): hl.or_else(mt['pc' + str(k + 1)], related_scores[mt.s].scores[k])
                                 for k in range(args.pc_num)})

        plot_pcs(mt, args)

        mt = mt.checkpoint(pcs_path, overwrite=True)
    else:
        logging.info("Detected matrix table with principal components annotated exists, loading that.")
        mt = hl.read_matrix_table(pcs_path)

    return mt


def annotate_with_external_data(mt, args):
    """Annotates variants with external datasets like CADD and gnomAD."""
    annotated_path = os.path.join(args.out_dir, "annotated.mt")

    if not utils.check_exists(annotated_path) or args.force:
        if args.mpc_ht:
            logging.info("Annotating variants with MPC info.")
            mt = va.annotate_variants_mpc(mt, args.mpc_ht)
        if args.cadd_ht:
            logging.info("Annotating variants with CADD info.")
            mt = va.annotate_variants_cadd(mt, args.cadd_ht)
        if args.gnomad_ht:
            logging.info("Annotating variants with Gnomad info.")
            mt = va.annotate_variants_gnomad(mt, args.gnomad_ht)
        if args.gnomad_mismatch_ht:
            logging.info("Annotating variants with gnomad mismatch info.")
            mt = va.annotate_variants_gnomad_mismatch(mt, args.gnomad_mismatch_ht)

        mt = mt.checkpoint(annotated_path, overwrite=True)
    else:
        mt = hl.read_matrix_table(annotated_path)

    return mt


def export_rows_cols(mt, args):
    """Export rows and columns of matrix table to text files."""
    mtcols = mt.cols()
    mtcols = mtcols.flatten()
    mtcols.export(os.path.join(args.out_dir, args.out_name + '_final_dataset_cols.tsv.gz'))

    mtrows = mt.rows()
    mtrows = mtrows.flatten()
    mtrows = mtrows.key_by().drop("vep.input")
    mtrows.export(os.path.join(args.out_dir, args.out_name + '_final_dataset_rows.tsv.gz'))


def main():
    hl.init()

    args = parse_arguments(sys.argv[1:])
    initialize_logger(args)
    verify_inputs(args)

    mt = load_and_annotate_samples(args)
    mt = run_variant_qc(mt, args, "low_pass")
    plot_variant_stats(mt, "low_pass")

    mt, mt_downsampled = calculate_relatedness(mt, args, "low_pass")
    mt, mt_downsampled = find_population_outliers(mt, mt_downsampled, args)
    mt = annotate_variants(mt, args)
    mt = impute_sex(mt, args, "low_pass")

    mt = run_samples_qc(mt, args)
    mt = run_variant_qc(mt, args, "final")

    mt = calculate_final_pcs(mt, args, "final")
    mt = annotate_with_external_data(mt, args)

    mt.write(os.path.join(args.out_dir, "final_qc.mt"), overwrite=True)
    utils.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.plot_folder)


if __name__ == "__main__":
    main()
