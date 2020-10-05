"""
Script to find and report de novo mutations, given parental relationship information.
"""
import argparse
import os
import time
import logging
import sys
from collections import OrderedDict
import utils


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
    exists = utils.check_if_gcloud_object_exists(kinship_annotated_fam)

    if (exists == 0) and (args.force == False):
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
        exists = utils.check_if_gcloud_object_exists(new_fam)

        if exists == 0:
            logging.info(f"Detected validated fam file exists, loading that: {new_fam}")
            new_fam = hl.import_table(new_fam)
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

    ############################################
    # Annotate dataset with gnomad frequencies #
    ############################################
    mt = va.annotate_variants_gnomad(mt, args.gnomad_ht)

    # First of all, use prior AF if it exists in Gnomad
    logging.info(f"Population for population allele frequency prior, from Gnomad: {args.gnomad_population}")
    mt = mt.annotate_rows(denovo_prior=hl.cond(
        (hl.len(mt.gnomad_filters) == 0) &
        hl.is_defined(mt.gnomad_freq[mt.gnomad_freq_index_dict[args.gnomad_population]]),
        mt.gnomad_freq[mt.gnomad_freq_index_dict[args.gnomad_population]].AF, 0))

    check = mt.aggregate_rows(hl.agg.counter(hl.is_defined(mt.denovo_prior)))
    logging.info(f"Missing AFs after adding gnomad AFs: {check}")

    prior_stats = mt.aggregate_rows(hl.agg.stats(mt.denovo_prior))
    logging.info(f"De novo prior AF stats: {prior_stats}")

    mt = mt.checkpoint(args.output_stem + "_gnomad_annotated_tmp.mt")

    #####################################################
    # Get de novo mutations, annotate with variant info #
    #####################################################
    denovos = hl.de_novo(mt, pedigree, pop_frequency_prior=mt.denovo_prior)
    denovos = denovos.checkpoint(args.output_stem + "_denovos_unannotated_tmp.ht")

    h.remove_preemptibles(args.cluster_name)
    denovos = denovos.key_by(denovos.locus, denovos.alleles)
    denovos = denovos.checkpoint(args.output_stem + "_rekeyed_tmp.ht")
    h.add_preemptibles(args.cluster_name, args.num_2nd_workers)

    mtrows = mt.rows()

    denovos = denovos.annotate(
        final_failing_variant_qc=mtrows[denovos.locus, denovos.alleles].final_failing_variant_qc,
        final_no_failing_samples_varqc=mtrows[denovos.locus, denovos.alleles].final_no_failing_samples_varqc,
        gene=mtrows[denovos.locus, denovos.alleles].gene,
        LOF=mtrows[denovos.locus, denovos.alleles].LOF,
        missense=mtrows[denovos.locus, denovos.alleles].missense
    )

    ###################################
    # Summarize # of de novo variants #
    ###################################
    denovo_count = denovos.count()
    high_confidence = denovos.filter((denovos.confidence == "HIGH") & (denovos.prior < 0.01))

    logging.info(f"Number of de novo mutations, total: {denovo_count}")
    logging.info(f"Count of de high-confidence novo mutations (with AF < 1%): {high_confidence}")

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
    parser.add_argument("--output_name", required=True, type=str, help="Output name for files.")
    parser.add_argument("--scripts_dir", required=True, type=str, help="Directory containing scripts for this library.")
    parser.add_argument("--output_dir", required=True, type=str, help="Output directory for output files.")
    parser.add_argument("--log_dir", required=True, type=str, help="Output directory for logs.")
    parser.add_argument("--cluster_name", required=True, type=str, help="Name of cluster for adding secondary workers")
    parser.add_argument("--num_2nd_workers", type=int, default=20, help="Number of secondary workers to add")
    args = parser.parse_args()

    scripts = ["variant_annotation.py", "helper_scripts.py"]
    for script in scripts:
        hl.spark_context().addPyFile(os.path.join(args.scripts_dir, script))

    import variant_annotation as va
    import helper_scripts as h

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
    denovo_table.export(args.output_stem + "_denovo_variants.txt")

    ########################
    # Copy logs and finish #
    ########################
    logging.info("Denovo pipeline ran successfully. Copying logs and exiting.")
    h.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.log_dir)
