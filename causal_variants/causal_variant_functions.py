"""
Functions used for finding putative causal variants
Author: Lea Urpa
August 2021
"""
import utils
import logging
import hail as hl


def count_case_control_carriers(mt, checkpoint_name, pheno_col='is_case_final', female_col='is_female_imputed'):
    """
    Annotate het and hom var carriers for variants in matrix table.

    :param mt: matrix table to annotate/filter
    :param args: arguments for pheno column, etc
    :return: returns annotated matrix table
    """
    ############################################################
    # Get count of samples that are cases and controls, report #
    ############################################################
    logging.info("Annotating het and hom var carrier counts in controls to variants.")

    case_count = mt.aggregate_cols(hl.agg.count_where(mt[pheno_col] == True))
    control_count = mt.aggregate_cols(hl.agg.count_where(mt[pheno_col] == False))
    missing = mt.aggregate_cols(hl.agg.count_where(~hl.is_defined(mt[pheno_col])))

    logging.info(f"Number of controls in dataset: {control_count}")
    logging.info(f"Number of cases in dataset: {case_count}")
    logging.info(f"Samples missing case/control information: {missing}")
    if missing > 0:
        logging.info(f"Warning- samples missing case/control status will be generally ignored in this pipeline.")

    ##################################################
    # Annotate control/case het count + homvar count #
    ##################################################
    mt = mt.annotate_rows(control_het_count=
                          hl.agg.filter((mt[pheno_col] == False) & hl.is_defined(mt[pheno_col]),
                                        hl.agg.count_where(mt.GT.is_het())))
    mt = mt.annotate_rows(control_homvar_count=
                          hl.agg.filter((mt[pheno_col] == False) & hl.is_defined(mt[pheno_col]),
                                        hl.agg.count_where(mt.GT.is_hom_var())))
    mt = mt.annotate_rows(case_het_count=
                          hl.agg.filter((mt[pheno_col] == True) & hl.is_defined(mt[pheno_col]),
                                        hl.agg.count_where(mt.GT.is_het())))
    mt = mt.annotate_rows(case_homvar_count=
                          hl.agg.filter((mt[pheno_col] == True) & hl.is_defined(mt[pheno_col]),
                                        hl.agg.count_where(mt.GT.is_hom_var())))

    ################################################################################
    # Count control homvar females and het males for hemizygous mutations on X chr #
    ################################################################################
    mt = mt.annotate_rows(control_homvar_count_female=
                          hl.agg.filter((mt[pheno_col] == False) & hl.is_defined(mt[pheno_col]) &
                                        (mt[female_col] == True) & hl.is_defined(mt[female_col]),
                                        hl.agg.count_where(mt.GT.is_hom_var())))

    mt = mt.annotate_rows(control_het_count_male=
                          hl.agg.filter((mt[pheno_col] == False) & hl.is_defined(mt[pheno_col]) &
                                        (mt[female_col] == False) & hl.is_defined(mt[female_col]),
                                        hl.agg.count_where(mt.GT.is_het())))

    mt = mt.annotate_rows(case_homvar_count_female=
                          hl.agg.filter((mt[pheno_col] == True) & hl.is_defined(mt[pheno_col]) &
                                        (mt[female_col] == True) & hl.is_defined(mt[female_col]),
                                        hl.agg.count_where(mt.GT.is_hom_var())))

    mt = mt.annotate_rows(case_het_count_male=
                          hl.agg.filter((mt[pheno_col] == True) & hl.is_defined(mt[pheno_col]) &
                                        (mt[female_col] == False) & hl.is_defined(mt[female_col]),
                                        hl.agg.count_where(mt.GT.is_het())))

    mt = mt.checkpoint(checkpoint_name, overwrite=True)

    return mt


def annotate_variants(mt, args):
    """
    Annotates matrix table variants with CADD, MPC and gnomad info.
    :param mt: matrix table to annotate
    :param args: arguments for cadd, mpc, and gnomad hail table locations
    :return: returns annotated matrix table
    """
    # TODO change output to depend on initial input name, but in output directory
    temp_filename = args.output_stem + "_cadd_mpc_gnomad_annotated_tmp.mt"

    if utils.check_exists(temp_filename) and (args.force == False):
        logging.info(f"Detected that matrix table annotated with CADD, MPC and Gnomad exist: {temp_filename}. "
                     f"Loading this.")
        mt = hl.read_matrix_table(temp_filename)
    else:
        args.force = True
        logging.info("Annotating matrix table with CADD, MPC and Gnomad.")

        h.add_preemptibles(args.cluster_name, args.num_preemptibles)

        #####################################
        # Annotate variants with CADD + MPC #
        #####################################
        mt = va.annotate_variants_cadd(mt, args.cadd_ht)
        mt = va.annotate_variants_mpc(mt, args.mpc_ht)

        #######################################################
        # Annotate variants with Gnomad population + mismatch #
        #######################################################
        mt = va.annotate_variants_gnomad(mt, args.gnomad_ht)
        args.gnomad_idx = mt.gnomad_popmax_index_dict.take(1)[0][args.gnomad_population]
        mt = va.annotate_variants_gnomad_mismatch(mt, args.gnomad_mismatch_ht)

        mt = mt.checkpoint(temp_filename, overwrite=True)

        h.remove_preemptibles(args.cluster_name)

    return mt

def annotate_population_thresholds(mt, args):
    """
    Annotate with True/False categories based on population parameters for recessive and dominant inheritance models
    in Gnomad population counts and allele frequencies and in control counts in dataset.

    :param mt: matrix table to annotate
    :param args: arguments giving max allowed carrier counts, max allele frequencies, and gnomad population index
    :return: returns annotated matrix table
    """
    # TODO change output to depend on initial input name, but in output directory
    temp_filename = args.output_stem + "_gnomad_control_thresholds_annotated_tmp.mt"

    if utils.check_exists(temp_filename) and (args.force == False):
        logging.info(f"Detected file with boolean columns for variant fulfulling population criteria: {temp_filename}. "
                     f"Loading this file.")
        mt = hl.read_matrix_table(temp_filename)
    else:
        args.force = True
        logging.info("Annotating with boolean columns for whether variant fulfills population threshold criteria.")

        ############################
        # Pull rows and checkpoint #
        ############################
        rows = mt.rows()
        rows = rows.checkpoint(args.output_stem + "_rows_tmp.ht", overwrite=True)

        ####################################################
        # Annotate whether variant dominant rare in gnomad #
        ####################################################
        rows = rows.annotate(dominant_rare_gnomad=hl.cond(
            # If gnomad data is good,
            (hl.len(rows.gnomad_filters) == 0) & hl.is_defined(rows.gnomad_popmax[args.gnomad_idx]) &
            (rows.gnomad_mismatch == False),
            # Ask whether the gnomad AC is small and there are no homozygotes
            hl.cond(
                (rows.gnomad_popmax[args.gnomad_idx].AC <= args.max_allowed_carrier_dominant) &
                (rows.gnomad_popmax[args.gnomad_idx].homozygote_count == 0),
                True, False),
            # else return None
            hl.null(hl.tbool)))

        #####################################################
        # Annotate whether variant recessive rare in gnomad #
        #####################################################
        rows = rows.annotate(recessive_rare_gnomad=hl.cond(
            (hl.len(rows.gnomad_filters) == 0) & hl.is_defined(rows.gnomad_popmax[args.gnomad_idx]) &
            (rows.gnomad_mismatch == False),
            hl.cond(
                (rows.gnomad_popmax[args.gnomad_idx].homozygote_count <= args.max_allowed_homozygotes_recessive)
                & (rows.gnomad_popmax[args.gnomad_idx].AF <= args.gnomad_AF_cutoff_recessive),
                True, False),
            hl.null(hl.tbool)))

        ###################################################################################################################
        # Find which population is the max for each variant, find index of males and females for that population in freqs #
        ###################################################################################################################
        rows = rows.annotate(popmax_population=rows.gnomad_popmax[args.gnomad_idx].pop)
        rows = rows.annotate(popmax_female_index=
                              rows.gnomad_freq_index_dict[args.gnomad_population + "_" + rows.popmax_population + "_female"])
        rows = rows.annotate(popmax_male_index=
                              rows.gnomad_freq_index_dict[args.gnomad_population + "_" + rows.popmax_population + "_male"])

        ######################################################
        # Annotate whether variant hemizygous rare in gnomad #
        ######################################################
        rows = rows.annotate(hemizygous_rare_gnomad=hl.cond(
            (hl.len(rows.gnomad_filters) == 0) & hl.is_defined(rows.gnomad_popmax[args.gnomad_idx]) &
            (rows.gnomad_mismatch == False) & (rows.locus.contig == 'chrX'),
            hl.cond(
                (rows.gnomad_freq[rows.popmax_male_index].AC <= args.max_allowed_carrier_dominant) &
                (rows.gnomad_freq[rows.popmax_female_index].homozygote_count <= args.max_allowed_homozygotes_recessive) &
                (rows.gnomad_popmax[args.gnomad_idx].AF <= args.gnomad_AF_cutoff_recessive),
                True, False),
            hl.null(hl.tbool)
        ))

        ######################################################################################
        # Annotate whether variant dominant/recessive/hemizygous rare in controls in dataset #
        ######################################################################################
        rows = rows.annotate(dominant_rare_controls=hl.cond(
            (rows.control_het_count <= args.max_allowed_carrier_dominant) &
            (rows.control_homvar_count == 0),
            True, False))

        rows = rows.annotate(recessive_rare_controls=hl.cond(
            (rows.control_homvar_count <= args.max_allowed_homozygotes_recessive),
            True, False))

        rows = rows.annotate(hemizygous_rare_controls=hl.cond(
            (rows.locus.contig == 'chrX') &
            (rows.control_homvar_count_female <= args.max_allowed_homozygotes_recessive) &
            (rows.control_het_count_male <= args.max_allowed_carrier_dominant),
            True, False))

        #####################################################
        # Annotate matrix table rows from 'rows' hail table #
        #####################################################
        mt = mt.annotate_rows(dominant_rare_gnomad=rows[mt.locus, mt.alleles].dominant_rare_gnomad,
                              recessive_rare_gnomad=rows[mt.locus, mt.alleles].recessive_rare_gnomad,
                              hemizygous_rare_gnomad=rows[mt.locus, mt.alleles].hemizygous_rare_gnomad,
                              dominant_rare_controls=rows[mt.locus, mt.alleles].dominant_rare_controls,
                              recessive_rare_controls=rows[mt.locus, mt.alleles].recessive_rare_controls,
                              hemizygous_rare_controls=rows[mt.locus, mt.alleles].hemizygous_rare_controls)

        mt = mt.checkpoint(temp_filename, overwrite=True)

    dominant_rare_gnomad = mt.aggregate_rows(hl.agg.counter(mt.dominant_rare_gnomad))
    recessive_rare_gnomad = mt.aggregate_rows(hl.agg.counter(mt.recessive_rare_gnomad))
    hemizygous_rare_gnomad = mt.aggregate_rows(hl.agg.counter(mt.hemizygous_rare_gnomad))
    dominant_rare_controls = mt.aggregate_rows(hl.agg.counter(mt.dominant_rare_controls))
    recessive_rare_controls = mt.aggregate_rows(hl.agg.counter(mt.recessive_rare_controls))
    hemizygous_rare_controls = mt.aggregate_rows(hl.agg.counter(mt.hemizygous_rare_controls))

    logging.info(f"Number of variants that are dominant rare in gnomad: {dominant_rare_gnomad}")
    logging.info(f"Number of variants that are recessive rare in gnomad: {recessive_rare_gnomad}")
    logging.info(f"Number of variants that are hemizygous rare in gnomad: {hemizygous_rare_gnomad}")

    logging.info(f"Number of variants that are dominant rare in controls: {dominant_rare_controls}")
    logging.info(f"Number of variants that are recessive rare in controls: {recessive_rare_controls}")
    logging.info(f"Number of variants that are hemizygous rare in controls: {hemizygous_rare_controls}")

    return mt