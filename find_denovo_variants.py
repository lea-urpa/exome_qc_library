"""
Script to find and report de novo mutations, given parental relationship information.
"""
import argparse
import os
import time
import logging
import sys
from collections import OrderedDict


def validate_pedigree(fam, kin, args):
    """
    Imports a fam file and kinship file (.kin king output), and looks for PO relationships as predicted in the data
    and validates that the kinship values support the relationship. If they don't, it removes those lines from the
    fam file and writes a new version for use in later steps.

    :param fam: fam file
    :param kin: kinship file output from King (.kin file).
    :param args: arguments for output directory
    :return: fam file, same as input if no errors, otherwise the new fam file with offending lines removed.
    """
    kinship_annotated_fam = args.output_stem + "_kinship_annotated_fam.txt"

    if utils.check_exists(kinship_annotated_fam) and (args.force == False):
        logging.info(f"Detected kinship annotated fam exists, loading file: {kinship_annotated_fam}")
        fam_ht = hl.import_table(kinship_annotated_fam, impute=True)
    else:
        args.force = True
        # Import pedigree file from fam as Hail table
        fam_ht = hl.import_table(fam, no_header=True)

        # Collect MID and PID as array, explode so one line with ID1 and FID or MID
        fam_ht = fam_ht.annotate(parent=[fam_ht.f2, fam_ht.f3])
        fam_ht = fam_ht.explode('parent')

        # Create key to merge with kinship data, explode (IDs may be switched in kinship data)
        fam_ht = fam_ht.annotate(kinship_key=[fam_ht.f0 + ":" + fam_ht.f1 + ":" + fam_ht.parent,
                                              fam_ht.f0 + ":" + fam_ht.parent + ":" + fam_ht.f1])
        fam_ht = fam_ht.explode('kinship_key')
        fam_ht = fam_ht.key_by('kinship_key')

        # Import King output files
        kinship = hl.import_table(kin, types={'Kinship': hl.tfloat64, 'IBS0': hl.tfloat64})
        kinship = kinship.annotate(kinship_key=kinship.FID + ":" + kinship.ID1 + ":" + kinship.ID2)
        kinship = kinship.key_by('kinship_key')

        # Add kinship information to fam ht
        fam_ht = fam_ht.annotate(Kinship=kinship.index(fam_ht.key).Kinship,
                                 IBS0=kinship.index(fam_ht.key).IBS0)

        # Export to tsv a file with the relationships and the kinship values side by side
        logging.info("Exporting kinship-annotated fam file")
        fam_ht.export(kinship_annotated_fam)

        # Remove rows in fam ht that do not have kinship data
        fam_ht = fam_ht.filter(hl.is_defined(fam_ht.Kinship))

        # Find PO relationships not supported by the data
        fam_ht = fam_ht.annotate(incorrect_po=hl.cond(
            (fam_ht.Kinship < 0.177) | (fam_ht.Kinship > 0.354) | (fam_ht.IBS0 > 0.1),
            True, False))

        failing_count = fam_ht.aggregate(hl.agg.count_where(fam_ht.incorrect_po == True))
        logging.info(f"Number of reported parent-offspring pairs unsupported by kinship or IBS0 values: {failing_count}")

        # If there are PO relationships not supported by the data, write new fam file excluding those individuals
        if failing_count > 0:
            # Print offending lines
            logging.info(fam_ht.filter(fam_ht.incorrect_po == True).show(failing_count))

            new_fam = fam.replace(".fam", "") + "_without_failing_PO_trios.fam"

            if utils.check_exists(new_fam):
                logging.info(f"Detected validated fam file exists, skipping writing new fam file.")
            else:
                # Create key of family ID, MID, PID
                fam_ht = fam_ht.annotate(fid_iid_mid_pid=fam_ht.f0 + ":" + fam_ht.f1 + ":" + fam_ht.f2 + ":" + fam_ht.f3)

                # Collect set of lines in fam ht that are failing kinship checks in either parent
                failing_fam_lines = fam_ht.aggregate(hl.agg.filter(fam_ht.incorrect_po == True,
                    hl.agg.collect_as_set(fam_ht.fid_iid_mid_pid)))

                # Filter out lines with incorrect PO from fam
                fam_ht = fam_ht.filter(hl.array(failing_fam_lines).contains(fam_ht.fid_iid_mid_pid), keep=False)

                # Drop extra columns and get unique lines
                fam_ht = fam_ht.key_by('f0', 'f1')
                fam_ht = fam_ht.drop('parent', 'kinship_key', 'Kinship', 'IBS0', 'incorrect_po', 'fid_iid_mid_pid')
                fam_ht = fam_ht.distinct()

                # Write to new fam file
                count = fam_ht.count()
                logging.info(f"Number of lines in new fam file, without lines containing PO relationships unsupported by kinship: "
                      f"{count}")

                fam_ht.export(new_fam)

            fam = new_fam

    return fam


def get_denovos(fam, mt, args):
    """
    Find denovos with Hail's hl.de_novo function

    :param fam: family pedigree to read in
    :param mt: matrix table with variants to find de novos
    :param args: argument for gnomad population to pull from population frequency dictionary of annotations
    :return: returns table of de novo results
    """
    h.add_preemptibles(args.cluster_name, args.num_2nd_workers)
    #########################################################
    # Read in trios, filter to just samples in matrix table #
    #########################################################
    pedigree = hl.Pedigree.read(fam)
    mt_samples = mt.s.take(mt.count_cols())
    pedigree = pedigree.filter_to(mt_samples)

    ####################################################
    # Report number of trios, number of complete trios #
    ####################################################
    trios = len(pedigree.trios)
    complete_trios = len(pedigree.complete_trios())
    logging.info(f"Number of trios: {trios}")
    logging.info(f"Number of complete trios: {complete_trios}")

    #############################
    # Pull base name and folder #
    #############################
    mt_name = args.mt
    if mt_name.endswith("/"):
        mt_name = mt_name.rstrip("/")

    folder, basename = os.path.split(mt_name)
    basename = basename.rstrip(".mt")

    ############################################
    # Annotate dataset with gnomad frequencies #
    ############################################
    if not utils.check_exists(os.path.join(folder, f"{basename}_gnomad_annotated.mt")):
        mt = va.annotate_variants_gnomad(mt, args.gnomad_ht)
        mt = mt.checkpoint(os.path.join(folder, f"{basename}_gnomad_annotated.mt"), overwrite=True)
    else:
        logging.info("Detected gnomad annotated matrix table exists. Loading that.")
        mt = hl.read_matrix_table(os.path.join(folder, f"{basename}_gnomad_annotated.mt"))

    # A little exposition
    #
    # Given a list of population priors, Hail uses the max of
    #  - n alt alleles - 1 / total alleles
    #  - given pop frequency
    #  - min pop prior (100/30000000)
    #
    # unless you give ignore_in_sample_allele_frequency, in which case it skips the first measure
    # For our data, which includes failing genotypes, the in-sample AF might be off
    # We will then filter failing genotypes, calculate n alt alleles -1 / total alleles
    # then annotate unfiltered dataset with those values
    # Then take the max of gnomad AFs and genotype-filtered sample AFs and set that as the gnomad prior
    # Then run de novo calling with the ignore_in_sample_allele_frequency flag.
    force = False

    # Filter genotypes
    if not utils.check_exists(os.path.join(args.output_dir, f"{args.output_name}_genotypes_filtered_tmp.mt")):
        mt_genofilt = mt.filter_entries((hl.len(mt.final_failing_depth_quality) == 0) &
                                        (hl.len(mt.final_failing_ab) == 0))
        mt_genofilt = mt_genofilt.checkpoint(os.path.join(args.output_dir, f"{args.output_name}_genotypes_filtered_tmp.mt"),
                                             overwrite=True)
        force = True
    else:
        mt_genofilt = hl.read_matrix_table(os.path.join(args.output_dir, f"{args.output_name}_genotypes_filtered_tmp.mt"))

    # Calculate in-sample allele frequency, and AF if gnomad N is included
    if (not utils.check_exists(os.path.join(args.output_dir, f"{args.output_name}_genotypes_filtered_rows.ht/"))) or \
            (force is True):
        # https://gnomad.broadinstitute.org/faq 56885 is gnomad v2 sample size
        gnomad_fin_AN = 2 * 56885
        mt_genofilt = mt_genofilt.annotate_rows(n_alt_alleles=hl.agg.sum(mt_genofilt.GT.n_alt_alleles()),
                                                total_alleles=2 * hl.agg.sum(hl.is_defined(mt_genofilt.GT)))

        mt_genofilt = mt_genofilt.annotate_rows(
            site_freq_filtered=(mt_genofilt.n_alt_alleles - 1) / mt_genofilt.total_alleles,
            site_freq_gnomad_n=(mt_genofilt.n_alt_alleles - 1) / (mt_genofilt.total_alleles + gnomad_fin_AN))

        genofilt_rows = mt_genofilt.rows()
        genofilt_rows = genofilt_rows.checkpoint(os.path.join(args.output_dir, f"{args.output_name}_genotypes_filtered_rows.ht/"),
                                                 overwrite=True)
        force = True
    else:
        genofilt_rows = hl.read_table(os.path.join(args.output_dir, f"{args.output_name}_genotypes_filtered_rows.ht/"))

    # Annotate original dataset with allele frequencies from genotype filtered dataset
    if (not utils.check_exists(os.path.join(args.output_dir, f"{args.output_name}_denovo_prior_annotated.mt/"))) or \
        (force is True):
        # Annotate alt alleles, allele frequencies from genotype filtered dataset to unfiltered datset
        mt = mt.annotate_rows(n_alt_alleles=genofilt_rows[mt.locus, mt.alleles].n_alt_alleles,
                              total_alleles=genofilt_rows[mt.locus, mt.alleles].total_alleles,
                              site_freq_filtered=genofilt_rows[mt.locus, mt.alleles].site_freq_filtered,
                              site_freq_gnomad_n=genofilt_rows[mt.locus, mt.alleles].site_freq_gnomad_n)

        # Pull out gnomad fin allele frequencies
        gnomad_population = "gnomad_fin"
        mt = mt.annotate_rows(gnomad_af=hl.if_else((hl.len(mt.gnomad_filters) == 0) & hl.is_defined(mt.gnomad_freq),
                                                   mt.gnomad_freq[mt.gnomad_freq_index_dict[gnomad_population]].AF,
                                                   0))
        mt = mt.annotate_rows(gnomad_af=hl.or_else(mt.gnomad_af, 0))

        check = mt.aggregate_rows(hl.agg.counter(hl.is_defined(mt.gnomad_af)))
        logging.info(f"Sanity check: gnomad_af is_defined count: {check}")

        # Annotate de novo prior as max of gnomad AF and site freq, with gnomad AN in denominator
        mt = mt.annotate_rows(denovo_prior=hl.max(mt.gnomad_af, mt.site_freq_gnomad_n))

        mt = mt.checkpoint(os.path.join(args.output_dir, f"{args.output_name}_denovo_prior_annotated.mt/"),
                           overwrite=True)
        force = True
    else:
        logging.info('Detected denovo prior annotated file exists, loading that')
        mt = hl.read_matrix_table(os.path.join(args.output_dir, f"{args.output_name}_denovo_prior_annotated.mt/"))

    #####################################################
    # Get de novo mutations, annotate with variant info #
    #####################################################
    # Pull rows to annotate de novos table
    if (not utils.check_exists(os.path.join(args.output_dir, f"{args.output_name}_rows.ht"))) or (force is True):
        rows = mt.rows()
        rows = rows.checkpoint(os.path.join(args.output_dir, f"{args.output_name}_rows.ht"), overwrite=True)
        force = True
    else:
        rows = hl.read_table(os.path.join(args.output_dir, f"{args.output_name}_rows.ht"))

    # Call de novos
    if (not utils.check_exists(os.path.join(args.output_dir, f"{args.output_name}_denovos_unannotated_tmp.ht"))) or \
            (force is True):
        denovos = hl.de_novo(mt, pedigree, pop_frequency_prior=mt.denovo_prior,
                             ignore_in_sample_allele_frequency=True)
        denovos = denovos.checkpoint(os.path.join(args.output_dir, f"{args.output_name}_denovos_unannotated_tmp.ht"),
                                     overwrite=True)
        force = True
    else:
        denovos = hl.read_table(os.path.join(args.output_dir, f"{args.output_name}_denovos_unannotated_tmp.ht"))

    #######################
    ## Annotate de novos ##
    #######################
    # Annotate de novos with variant information
    if (not utils.check_exists(os.path.join(args.output_dir, f"{args.output_name}_denovos_annotated_tmp.ht"))) or \
            (force is True):
        denovos = denovos.key_by(denovos.locus, denovos.alleles)

        denovos = denovos.annotate(
            final_failing_variant_qc=rows[denovos.locus, denovos.alleles].final_failing_variant_qc,
            final_no_failing_samples_varqc=rows[denovos.locus, denovos.alleles].final_no_failing_samples_varqc,
            gene=rows[denovos.locus, denovos.alleles].gene,
            gene_most_severe_conseq=rows[denovos.locus, denovos.alleles].gene_most_severe_conseq,
            LOF=rows[denovos.locus, denovos.alleles].LOF,
            missense=rows[denovos.locus, denovos.alleles].missense,
            synonymous=rows[denovos.locus, denovos.alleles].synonymous,
            vep=rows[denovos.locus, denovos.alleles].vep)
        force = True
    else:
        denovos = hl.read_table(os.path.join(args.output_dir, f"{args.output_name}_denovos_annotated_tmp.ht"))

    # Annotate with CADD and MPC, if given
    if (not utils.check_exists(os.path.join(args.output_dir, f"{args.output_name}_denovos_annotated_tmp2.ht"))) or \
            (force is True):
        # Annotate CADD + MPC values
        if args.cadd_ht is not None:
            CADD = hl.read_table(args.cadd_ht)
            denovos = denovos.annotate(CADD_phred=CADD[denovos.locus, denovos.alleles].PHRED)

        if args.mpc_ht is not None:
            MPC = hl.read_table(args.mpc_ht)
            denovos = denovos.annotate(MPC=MPC[denovos.locus, denovos.alleles].MPC)
            denovos = denovos.checkpoint(os.path.join(args.output_dir, f"{args.output_name}_denovos_annotated_tmp2.ht"),
                                         overwrite=True)

        if args.segdup_intervals_file is not None:
            segdup_intervals = hl.import_locus_intervals(args.segdup_intervals_file, reference_genome=args.reference_genome,
                                                         skip_invalid_intervals=True)
            denovos = denovos.annotate(in_segdup_region=hl.if_else(hl.is_defined(segdup_intervals[denovos.locus]),
                                                                   True, False))

        force = True
    else:
        denovos = hl.read_table(os.path.join(args.output_dir, f"{args.output_name}_denovos_annotated_tmp2.ht"))

    # Annotate genes with gene list and constraint values

    if args.gnomad_gene_metrics is not None:
        if (not utils.check_exists(os.path.join(args.output_dir, f"{args.output_name}_denovos_annotated_tmp3.ht"))) or \
                (force is True):
            # Pull de novo genes and explode by gene
            denovo_genes = denovos.select(denovos.id, denovos.gene)
            denovo_genes = denovo_genes.explode(denovo_genes.gene)

            # Annotate pLI metrics for genes
            gene_metrics = hl.import_table(args.gnomad_gene_metrics, types={'pLI': hl.tfloat64}, key='gene')
            denovo_genes = denovo_genes.annotate(pLI=gene_metrics[denovo_genes.gene].pLI)

            vars_missing_pLI = denovo_genes.aggregate(hl.agg.counter(hl.is_defined(denovo_genes.pLI)))
            print(f"Count of variants where pLI values are missing (False) or not (True): {vars_missing_pLI}")

            gene_count = denovo_genes.group_by("gene").aggregate(pLI=hl.agg.mean(denovo_genes.pLI))
            missing_pLI = gene_count.aggregate(hl.agg.counter(hl.is_defined(gene_count.pLI)))
            print(f"Count of genes where pLI values are missing (False) or not (True): {missing_pLI}")

            # Group genes, pLI metrics and annotate to main de novos table
            denovo_vars = denovo_genes.group_by('locus', 'alleles', 'id').aggregate(
                pLI=hl.array(hl.agg.collect_as_set(denovo_genes.pLI)))
            denovo_vars = denovo_vars.key_by('locus', 'alleles', 'id')

            denovos = denovos.key_by('locus', 'alleles', 'id')
            denovos = denovos.annotate(pLI=denovo_vars[denovos.locus, denovos.alleles, denovos.id].pLI)

            denovos = denovos.checkpoint(os.path.join(args.output_dir, f"{args.output_name}_denovos_annotated_tmp3.ht"),
                                         overwrite=True)
        else:
            denovos = hl.read_table(os.path.join(args.output_dir, f"{args.output_name}_denovos_annotated_tmp3.ht"))

    # maybe rekeying fails above? use this in that case
    #h.remove_preemptibles(args.cluster_name)
    #h.add_preemptibles(args.cluster_name, args.num_2nd_workers)

    return denovos


if __name__ == "__main__":
    import hail as hl
    hl.init()

    ###################
    # Parse arguments #
    ###################
    parser = argparse.ArgumentParser(description="v9 exome sequencing dataset quality control pipeline.")
    parser.add_argument("-mt", type=str, help="Input matrix table to run analysis on")
    parser.add_argument("--fam", type=str, required=True, help="Pedigree file in Plink .fam format.")
    parser.add_argument("--kin", type=str, required=True,
                        help="Kinship output from King, .kin file, corresponding to fam.")
    parser.add_argument("--gnomad_ht", required=True, type=str, help="File name of gnomad hail table.")
    parser.add_argument("--gnomad_population", type=str, default="gnomad_fin",
                        choices=[])
    parser.add_argument("--reference_genome", default="GRCh38")
    parser.add_argument("--cadd_ht", type=str, help="Location of CADD hail table")
    parser.add_argument("--mpc_ht", type=str, help="Location of MPC hail table")
    parser.add_argument("--segdup_intervals_file", type=str, help="Location of segmental duplication intervals file.")
    parser.add_argument("--gnomad_gene_metrics", type=str, help="Location of gnomad gene constraint metrics, .bgz file")
    parser.add_argument("--output_name", required=True, type=str, help="Output name for files.")
    parser.add_argument("--scripts_dir", required=True, type=str, help="Directory containing scripts for this library.")
    parser.add_argument("--output_dir", required=True, type=str, help="Output directory for output files.")
    parser.add_argument("--log_dir", required=True, type=str, help="Output directory for logs.")
    parser.add_argument("--cluster_name", required=True, type=str, help="Name of cluster for adding secondary workers")
    parser.add_argument("--num_2nd_workers", type=int, default=20, help="Number of secondary workers to add")
    args = parser.parse_args()

    scripts = ["variant_annotation.py", "helper_scripts.py", "utils.py"]
    for script in scripts:
        hl.spark_context().addPyFile(os.path.join(args.scripts_dir, script))

    import variant_annotation as va
    import helper_scripts as h
    import utils

    args.output_stem = os.path.join(args.output_dir, args.output_name)
    args.force = False

    ####################
    # Configure logger #
    ####################
    datestr = time.strftime("%Y.%m.%d")  # Used for output folder
    timestr = time.strftime("%Y.%m.%d-%H.%M.%S")  # Used for output files, for more than one run per day
    args.log_file = 'find-denovo-variants_' + timestr + '.txt'

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

    ########################
    # Read in matrix table #
    ########################
    qcd_mt = hl.read_matrix_table(args.mt)

    ###########################################
    # Validate parent/offspring relationships #
    ###########################################
    validated_fam = validate_pedigree(args.fam, args.kin, args)

    ########################
    # Find denovo variants #
    ########################
    denovo_table = get_denovos(validated_fam, qcd_mt, args)

    denovo_table.write(args.output_stem + "_denovo_variants.ht", overwrite=True)
    denovo_table = denovo_table.flatten()
    denovo_table = denovo_table.drop("vep.input")
    denovo_table.export(args.output_stem + "_denovo_variants.txt")

    ########################
    # Copy logs and finish #
    ########################
    logging.info("Denovo pipeline ran successfully. Copying logs and exiting.")
    h.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.log_dir)
