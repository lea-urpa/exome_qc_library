"""
Functions used for finding putative causal variants
Author: Lea Urpa
August 2021
"""
import utils
import logging
import hail as hl
import variant_annotation as va


def count_case_control_carriers(mt, checkpoint_name, pheno_col='is_case_final', female_col='is_female_imputed'):
    """
    Annotate het and hom var carriers for variants in matrix table.

    :param mt: matrix table to annotate
    :param checkpoint_name: checkpoint name
    :param pheno_col: column ID name for case control status (must be boolean)
    :param female_col: column ID name for female/male status (must be boolean, female = True)
    :return:
    """
    logging.info("Annotating het and hom var carrier counts in controls to variants.")
    annotations_to_transfer = [
        'control_het_count', 'control_homvar_count', 'case_het_count', 'case_homvar_count',
        'control_homvar_count_female', 'control_het_count_male', 'case_homvar_count_female',
        'case_het_count_male'
    ]

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

    return mt, annotations_to_transfer


def check_variants_annotated(mt, checkpoint_name, cadd_ht, mpc_ht, gnomad_ht, gnomad_mismatch_ht):
    """
    Checks if matrix table is already annotated with CADD, MPC, and gnomad, and annotates if not.

    :param mt: matrix table to annotate
    :param checkpoint_name: checkpoint
    :param cadd_ht:
    :param mpc_ht:
    :param gnomad_ht:
    :param gnomad_mismatch_ht:
    :return:
    """
    """
    Annotates matrix table variants with CADD, MPC and gnomad info.
    :param mt: matrix table to annotate
    :param args: arguments for cadd, mpc, and gnomad hail table locations
    :return: returns annotated matrix table
    """
    try:
        test = hl.is_defined(mt.CADD_PHRED)
        cadd_not_annot = False
    except:
        cadd_not_annot = True

    try:
        test = hl.is_defined(mt.MPC)
        mpc_not_annot = False
    except:
        mpc_not_annot = True

    try:
        test = hl.is_defined(mt.gnomad_freq)
        gnomad_not_annot = False
    except:
        gnomad_not_annot = True

    try:
        test = hl.is_defined(mt.gnomad_mismatch_pvalue)
        gnomad_mismatch_not_annot = False
    except:
        gnomad_mismatch_not_annot = True

    if cadd_not_annot:
        logging.info("Annotating matrix table with CADD")
        mt = va.annotate_variants_cadd(mt, cadd_ht)
    if mpc_not_annot:
        logging.info("Annotating matrix table with MPC")
        mt = va.annotate_variants_mpc(mt, mpc_ht)
    if gnomad_not_annot:
        logging.info("Annotating matrix table with gnomad variant information")
        mt = va.annotate_variants_gnomad(mt, gnomad_ht)
    if gnomad_mismatch_not_annot:
        mt = va.annotate_variants_gnomad_mismatch(mt, gnomad_mismatch_ht)

    if any(cadd_not_annot, mpc_not_annot, gnomad_not_annot, gnomad_mismatch_not_annot):
        mt = mt.checkpoint(checkpoint_name, overwrite=True)

    return mt


def annotate_control_thresholds(mt, checkpoint_name, gnomad_population, max_allowed_carrier_dominant=10,
                                   max_allowed_homozygotes_recessive=10, gnomad_AF_cutoff_recessive=0.01,
                                   ):
    """
    Annotate with True/False categories based on population parameters for recessive and dominant inheritance models
    in Gnomad population counts and allele frequencies and in control counts in dataset.

    :param mt: matrix table to annotate
    :param args: arguments giving max allowed carrier counts, max allele frequencies, and gnomad population index
    :return: returns annotated matrix table
    """
    logging.info("Annotating with boolean columns for whether variant fulfills population threshold criteria.")

    ############################
    # Pull rows and checkpoint #
    ############################
    rows = mt.rows()
    rows = rows.checkpoint(checkpoint_name + "_rows_tmp.ht", overwrite=True)

    gnomad_idx = mt.gnomad_popmax_index_dict.take(1)[0][gnomad_population]

    ####################################################
    # Annotate whether variant dominant rare in gnomad #
    ####################################################
    rows = rows.annotate(dominant_rare_gnomad=hl.cond(
        # If gnomad data is good,
        (hl.len(rows.gnomad_filters) == 0) & hl.is_defined(rows.gnomad_popmax[gnomad_idx]) &
        (rows.gnomad_mismatch == False),
        # Ask whether the gnomad AC is small and there are no homozygotes
        hl.cond(
            (rows.gnomad_popmax[gnomad_idx].AC <= max_allowed_carrier_dominant) &
            (rows.gnomad_popmax[gnomad_idx].homozygote_count == 0),
            True, False),
        # else return None
        hl.null(hl.tbool)))

    #####################################################
    # Annotate whether variant recessive rare in gnomad #
    #####################################################
    rows = rows.annotate(recessive_rare_gnomad=hl.cond(
        (hl.len(rows.gnomad_filters) == 0) & hl.is_defined(rows.gnomad_popmax[gnomad_idx]) &
        (rows.gnomad_mismatch == False),
        hl.cond(
            (rows.gnomad_popmax[gnomad_idx].homozygote_count <= max_allowed_homozygotes_recessive)
            & (rows.gnomad_popmax[gnomad_idx].AF <= gnomad_AF_cutoff_recessive),
            True, False),
        hl.null(hl.tbool)))

    ###################################################################################################################
    # Find which population is the max for each variant, find index of males and females for that population in freqs #
    ###################################################################################################################
    rows = rows.annotate(popmax_population=rows.gnomad_popmax[gnomad_idx].pop)
    rows = rows.annotate(popmax_female_index=
                          rows.gnomad_freq_index_dict[gnomad_population + "_" + rows.popmax_population + "_female"])
    rows = rows.annotate(popmax_male_index=
                          rows.gnomad_freq_index_dict[gnomad_population + "_" + rows.popmax_population + "_male"])

    ######################################################
    # Annotate whether variant hemizygous rare in gnomad #
    ######################################################
    rows = rows.annotate(hemizygous_rare_gnomad=hl.cond(
        (hl.len(rows.gnomad_filters) == 0) & hl.is_defined(rows.gnomad_popmax[gnomad_idx]) &
        (rows.gnomad_mismatch == False) & (rows.locus.contig == 'chrX'),
        hl.cond(
            (rows.gnomad_freq[rows.popmax_male_index].AC <= max_allowed_carrier_dominant) &
            (rows.gnomad_freq[rows.popmax_female_index].homozygote_count <= max_allowed_homozygotes_recessive) &
            (rows.gnomad_popmax[gnomad_idx].AF <= gnomad_AF_cutoff_recessive),
            True, False),
        hl.null(hl.tbool)
    ))

    ######################################################################################
    # Annotate whether variant dominant/recessive/hemizygous rare in controls in dataset #
    ######################################################################################
    rows = rows.annotate(dominant_rare_controls=hl.cond(
        (rows.control_het_count <= max_allowed_carrier_dominant) &
        (rows.control_homvar_count == 0),
        True, False))

    rows = rows.annotate(recessive_rare_controls=hl.cond(
        (rows.control_homvar_count <= max_allowed_homozygotes_recessive),
        True, False))

    rows = rows.annotate(hemizygous_rare_controls=hl.cond(
        (rows.locus.contig == 'chrX') &
        (rows.control_homvar_count_female <= max_allowed_homozygotes_recessive) &
        (rows.control_het_count_male <= max_allowed_carrier_dominant),
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

    mt = mt.checkpoint(checkpoint_name, overwrite=True)

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