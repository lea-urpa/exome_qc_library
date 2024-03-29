"""
Functions to annotate variants in a matrix table.

Created by Lea Urpa in 2019, some from Mitja's original code.
"""
import logging
import hail as hl


def annotate_variants(mt):
    '''
    Takes matrix table and annotates variants with gene, LOF and missense annotations by parsing VEP annotations.

    :param mt: matrix table to annotate
    :return: returns matrix table with new row annotations gene, LOF, and missense.
    '''
    try:
        test = hl.is_defined(mt.row.was_split)
    except Exception as e:
        print('Split multi-allelics before running!')
        print(e)
        return

    # If there is no canonical and protein-coding transcript consequence for that variant,
    # give the gene corresponding to the most severe consequence.
    # If there is a canonical and protein-coding transcript consequence for that variant,
    # give the gene symbol associated with that transcript consequence.
    # Also, if the consequence terms do not contain upstream or downstream gene variant, these can still be canonical.
    canon_pc = mt.row.vep.transcript_consequences.filter(
        lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding') &
                  (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                  ~x.consequence_terms.contains("downstream_gene_variant") &
                  ~x.consequence_terms.contains("upstream_gene_variant") &
                  ~x.consequence_terms.contains("intron_variant")
    )

    most_severe = mt.row.vep.transcript_consequences.filter(
        lambda x: x.consequence_terms.contains(mt.row.vep.most_severe_consequence) &
                  (x.biotype == "protein_coding") & (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                  ~x.consequence_terms.contains("downstream_gene_variant") &
                  ~x.consequence_terms.contains("upstream_gene_variant") &
                  ~x.consequence_terms.contains("intron_variant")
    )

    mt = mt.annotate_rows(
        gene=hl.if_else(hl.any(
            lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding') &
                      (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                      ~x.consequence_terms.contains("downstream_gene_variant") &
                      ~x.consequence_terms.contains("upstream_gene_variant") &
                      ~x.consequence_terms.contains("intron_variant"),
            mt.row.vep.transcript_consequences),
            canon_pc.map(lambda x:
                         hl.struct(gene_symbol=x.gene_symbol, gene_id=x.gene_id, consequence=x.consequence_terms)),
            most_severe.map(lambda x:
                            hl.struct(gene_symbol=x.gene_symbol, gene_id=x.gene_id, consequence=x.consequence_terms))),
    missing_false=True)

    # Same thing but for variants predicted to be noncoding in that transcript
    canon_pc = mt.row.vep.transcript_consequences.filter(
        lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding') &
                  (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                  (x.consequence_terms.contains("downstream_gene_variant") |
                   x.consequence_terms.contains("upstream_gene_variant") |
                   x.consequence_terms.contains("intron_variant"))
    )

    most_severe = mt.row.vep.transcript_consequences.filter(
        lambda x: x.consequence_terms.contains(mt.row.vep.most_severe_consequence) &
                  (x.biotype == "protein_coding") & (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                  (x.consequence_terms.contains("downstream_gene_variant") |
                   x.consequence_terms.contains("upstream_gene_variant") |
                   x.consequence_terms.contains("intron_variant"))
    )

    mt = mt.annotate_rows(
        gene_noncoding=hl.if_else(hl.any(
            lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding') &
                      (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                      (x.consequence_terms.contains("downstream_gene_variant") |
                       x.consequence_terms.contains("upstream_gene_variant") |
                       x.consequence_terms.contains("intron_variant")),
            mt.row.vep.transcript_consequences),
            canon_pc.map(lambda x:
                         hl.struct(gene_symbol=x.gene_symbol, gene_id=x.gene_id, consequence=x.consequence_terms)),
            most_severe.map(lambda x:
                            hl.struct(gene_symbol=x.gene_symbol, gene_id=x.gene_id, consequence=x.consequence_terms))),
    missing_false=True)

    # either if there is a canonical and protein coding transcript consequence for that variant,
    # and the lof annotation is not missing and equal to HC, and the lof flag is missing or is blank,
    # or if there isn't a canonical and protein coding transcript consequence for that variant and the
    # transcript consequence with consequence terms containing the most severe consequence term has lof not missing,
    # is equal to HC, and lof flags missing or blank,
    # true, else false

    canon_pc = mt.row.vep.transcript_consequences\
                         .filter(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding'))
    most_severe = mt.row.vep.transcript_consequences\
                            .filter(lambda x: x.consequence_terms.contains(
                                              mt.row.vep.most_severe_consequence))

    canon_bool = (hl.any(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding'),
                         mt.row.vep.transcript_consequences) &
                  hl.any(lambda x: hl.is_defined(x.lof), canon_pc) &
                  (canon_pc.map(lambda x: x.lof) == ["HC"]) &
                  (hl.all(lambda x: hl.is_missing(x.lof_flags) | (x.lof_flags == ""), canon_pc)))

    non_canon_bool = (~(hl.any(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding'),
                               mt.row.vep.transcript_consequences)) &
                      hl.any(lambda x: hl.is_defined(x.lof), most_severe) &
                      (most_severe.map(lambda x: x.lof) == ["HC"]) &
                      (hl.all(lambda x: hl.is_missing(x.lof_flags) | (x.lof_flags == ""), most_severe)))

    mt = mt.annotate_rows(LOF=hl.if_else(canon_bool | non_canon_bool, True, False))

    # Either if there is a canonical and protein coding transcript consequence for that variant
    # whose consequence terms contain "missense variant"
    # or if there is not a canonical and protein coding transcript consequence for that variant,
    # but the most severe consequence is "missense variant"
    # or if if there is a canonical and protein coding transcript consequence for that variant
    # whose consequence terms contain "inframe deletion"
    # or if there is not a canonical and protein coding transcript consequence for that variant,
    # but the variant's most severe consequence is "inframe deletion"
    # true else false

    canon_pc = mt.row.vep.transcript_consequences\
                         .filter(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding'))

    canon_missense_bool = canon_pc.map(lambda x: x.consequence_terms).contains(["missense_variant"])
    noncanon_missense_bool = (~(hl.any(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding'),
                                mt.row.vep.transcript_consequences)) &
                              (mt.row.vep.most_severe_consequence == "missense_variant"))

    canon_inframe_bool = canon_pc.map(lambda x: x.consequence_terms).contains(["inframe_deletion"])
    noncanon_inframe_bool = (~(hl.any(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding'),
                               mt.row.vep.transcript_consequences)) &
                             (mt.row.vep.most_severe_consequence == "inframe_deletion"))

    canon_inframe_ins_bool = canon_pc.map(lambda x: x.consequence_terms).contains(["inframe_insertion"])
    noncanon_inframe_ins_bool = (~(hl.any(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding'),
                               mt.row.vep.transcript_consequences)) &
                             (mt.row.vep.most_severe_consequence == "inframe_insertion"))

    mt = mt.annotate_rows(missense=hl.if_else((canon_missense_bool | noncanon_missense_bool | canon_inframe_bool |
                                               noncanon_inframe_bool | canon_inframe_ins_bool |
                                               noncanon_inframe_ins_bool),
                                              True, False))

    # If the most severe consequence is "synonymous_variant", true else false
    mt = mt.annotate_rows(synonymous=hl.if_else(mt.row.vep.most_severe_consequence == "synonymous_variant", True, False))

    # When there is a transcript consequence for that variant that is canonical,
    # protein coding, and lof = "HC", its lof flags
    # When there is not a transcript consequence for that variant that is canonical and protein coding,
    # but there is a transcript consequence whose consequence terms contains the most severe consequence
    # and its lof == HC, its lof flags
    # else blank

    canon_bool = hl.any(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding'),
                        mt.row.vep.transcript_consequences)
    canon_hc_bool = hl.any(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding') &
                                     (x.lof == 'HC'), mt.row.vep.transcript_consequences)
    canon_pc_hc = mt.row.vep.transcript_consequences.filter(lambda x:
                                                            (x.canonical == 1) & (x.biotype == 'protein_coding') &
                                                            (x.lof == "HC"))
    most_severe_bool = hl.any(lambda x: (x.consequence_terms.contains(mt.row.vep.most_severe_consequence)) &
                                        (x.lof == 'HC'),
                              mt.row.vep.transcript_consequences)
    most_severe_hc = mt.row.vep.transcript_consequences.filter(lambda x:
                                                               (x.consequence_terms.contains(
                                                                   mt.row.vep.most_severe_consequence)) &
                                                               (x.lof == "HC"))

    mt = mt.annotate_rows(LOF_flag=hl.case().when(canon_hc_bool,
                                                  canon_pc_hc.map(lambda x: x.lof_flags)
                                           ).when(~canon_bool & most_severe_bool,
                                                  most_severe_hc.map(lambda x: x.lof_flags)
                                           ).default([""]))

    return mt


def sex_aware_variant_annotations(mt, pheno_col=None, prefix=""):
    '''
    Creates sex-aware variant annotations for call rate, allele count, and allele number.

    :param mt: matrix table to annotate
    :param sex_col: string referring to column in the matrix table giving sex information
    :param male_tag: string or boolean in the column referring to males
    :param female_tag: string or boolean in the column referring to females
    :return: Returns matrix table with new row annotations male_hets, male_homvars, male_calls, female_hets,
    female_homvars, female_calls, sexaware_call_rate, sexaware_ac and sexaware_an.
    '''
    ####################################
    # Do sex-aware variant annotations #
    ####################################
    if (not prefix.endswith("_")) and (prefix != ""):
        prefix = prefix + "_"

    male_hets = prefix + "male_hets"
    male_homvars = prefix + "male_homvars"
    male_calls = prefix + "male_calls"
    female_hets = prefix + "female_hets"
    female_homvars = prefix + "female_homvars"
    female_calls = prefix + "female_calls"
    sexaware_call_rate = prefix + "sexaware_call_rate"
    sexaware_AC = prefix + "sexaware_AC"
    sexaware_AN = prefix + "sexaware_AN"

    logging.info('Annotating sex-aware variant annotations.')
    is_female = (mt.is_female_imputed == True) & hl.is_defined(mt.is_female_imputed)
    is_male = (mt.is_female_imputed == False) & hl.is_defined(mt.is_female_imputed)

    num_males = mt.aggregate_cols(hl.agg.count_where(is_male))
    num_females = mt.aggregate_cols(hl.agg.count_where(is_female))

    mt = mt.annotate_rows(
        **{male_hets: hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) & is_male),
           male_homvars: hl.agg.count_where(mt.GT.is_hom_var() & hl.is_defined(mt.GT) & is_male),
           male_calls: hl.agg.count_where(hl.is_defined(mt.GT) & is_male),
           female_hets: hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) & is_female),
           female_homvars: hl.agg.count_where(mt.GT.is_hom_var() & hl.is_defined(mt.GT) & is_female),
           female_calls: hl.agg.count_where(hl.is_defined(mt.GT) & is_female)}
    )

    mt = mt.annotate_rows(
        **{sexaware_call_rate: (hl.case()
                                  .when(mt.locus.in_y_nonpar() & hl.is_defined(mt.locus),
                                        (mt[male_calls] / num_males))
                                  .when(mt.locus.in_x_nonpar() & hl.is_defined(mt.locus),
                                        (mt[male_calls] + 2*mt[female_calls]) / (num_males + 2*num_females))
                                  .default((mt[male_calls] + mt[female_calls]) / (num_males + num_females))),
           sexaware_AC: (hl.case()  # MINOR allele count
                           .when(mt.locus.in_y_nonpar(), mt[male_homvars])
                           .when(mt.locus.in_x_nonpar(), mt[male_homvars] + mt[female_hets] + 2*mt[female_homvars])
                           .default(mt[male_hets] + 2*mt[male_homvars] + mt[female_hets] + 2*mt[female_homvars])),
           sexaware_AN: (hl.case()
                           .when(mt.locus.in_y_nonpar(), mt[male_calls])
                           .when(mt.locus.in_x_nonpar(), mt[male_calls] + 2*mt[female_calls])
                           .default(2*mt[male_calls] + 2*mt[female_calls]))}
    )

    annotations_to_transfer = [sexaware_call_rate, sexaware_AC, sexaware_AN]

    ##################################################
    # Do case-specific sex-aware variant annotations #
    ##################################################
    if pheno_col is not None:
        logging.info("Annotating sex-aware variant annotations, taking case/control status into account.")
        case_female = (mt.is_female_imputed == True) & (mt[pheno_col] == True)
        case_male = (mt.is_female_imputed == False) & (mt[pheno_col] == True)

        num_case_females = mt.aggregate_cols(hl.agg.count_where(case_female))
        num_case_males = mt.aggregate_cols(hl.agg.count_where(case_male))

        case_male_hets = prefix + "case_male_hets"
        case_male_homvars = prefix + "case_male_homvars"
        case_male_calls = prefix + "case_male_calls"
        case_female_hets = prefix + "case_female_hets"
        case_female_homvars = prefix + "case_female_homvars"
        case_female_calls = prefix + "case_female_calls"
        sexaware_case_call_rate = prefix + "sexaware_case_call_rate"
        sexaware_case_AC = prefix + "sexaware_case_AC"
        sexaware_case_AN = prefix + "sexaware_case_AN"

        mt = mt.annotate_rows(
            **{case_male_hets: hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) & case_male),
               case_male_homvars: hl.agg.count_where(mt.GT.is_hom_var() & hl.is_defined(mt.GT) & case_male),
               case_male_calls: hl.agg.count_where(hl.is_defined(mt.GT) & case_male),
               case_female_hets: hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) & case_female),
               case_female_homvars: hl.agg.count_where(mt.GT.is_hom_var() & hl.is_defined(mt.GT) & case_female),
               case_female_calls:hl.agg.count_where(hl.is_defined(mt.GT) & case_female)}
        )

        mt = mt.annotate_rows(**{
            sexaware_case_call_rate:
                              (hl.case()
                                  .when(mt.locus.in_y_nonpar() & hl.is_defined(mt.locus),
                                        (mt[case_male_calls] / num_case_males))
                                  .when(mt.locus.in_x_nonpar() & hl.is_defined(mt.locus),
                                        (mt[case_male_calls] + 2 * num_case_females) /
                                        (num_case_males + 2 * num_case_females))
                                  .default((mt[case_male_calls] + mt[case_female_calls]) /
                                           (num_case_males + num_case_females))),
            sexaware_case_AC:
                              (hl.case()  # MINOR allele count
                                 .when(mt.locus.in_y_nonpar(), mt[case_male_homvars])
                                 .when(mt.locus.in_x_nonpar(),
                                       mt[case_male_homvars] + mt[case_female_hets] + 2 * mt[case_female_homvars])
                                  .default(mt[case_male_hets] + 2 * mt[case_male_homvars] +
                                           mt[case_female_hets] + 2 * mt[case_female_homvars])),
            sexaware_case_AN:
                              (hl.case()
                                 .when(mt.locus.in_y_nonpar(), mt[case_male_calls])
                                 .when(mt.locus.in_x_nonpar(), mt[case_male_calls] + 2 * mt[case_female_calls])
                                 .default(2 * mt[case_male_calls] + 2 * mt[case_female_calls]))})

        annotations_to_transfer.extend(
            [case_male_hets, case_male_homvars, case_male_calls, case_female_hets, case_female_homvars,
             case_female_calls, sexaware_case_call_rate, sexaware_case_AC, sexaware_case_AN])

        cont_female = (mt.is_female_imputed == True) & (mt[pheno_col] == False)
        cont_male = (mt.is_female_imputed == False) & (mt[pheno_col] == False)

        num_cont_females = mt.aggregate_cols(hl.agg.count_where(cont_female))
        num_cont_males = mt.aggregate_cols(hl.agg.count_where(cont_male))

        cont_male_hets = prefix + "cont_male_hets"
        cont_male_homvars = prefix + "cont_male_homvars"
        cont_male_calls = prefix + "cont_male_calls"
        cont_female_hets = prefix + "cont_female_hets"
        cont_female_homvars = prefix + "cont_female_homvars"
        cont_female_calls = prefix + "cont_female_calls"
        sexaware_cont_call_rate = prefix + "sexaware_cont_call_rate"
        sexaware_cont_AC = prefix + "sexaware_cont_AC"
        sexaware_cont_AN = prefix + "sexaware_cont_AN"

        mt = mt.annotate_rows(
            **{cont_male_hets: hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) & cont_male),
               cont_male_homvars: hl.agg.count_where(mt.GT.is_hom_var() & hl.is_defined(mt.GT) & cont_male),
               cont_male_calls: hl.agg.count_where(hl.is_defined(mt.GT) & cont_male),
               cont_female_hets: hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) & cont_female),
               cont_female_homvars: hl.agg.count_where(mt.GT.is_hom_var() & hl.is_defined(mt.GT) & cont_female),
               cont_female_calls: hl.agg.count_where(hl.is_defined(mt.GT) & cont_female)})

        mt = mt.annotate_rows(
            **{sexaware_cont_call_rate:
                              (hl.case()
                               .when(mt.locus.in_y_nonpar() & hl.is_defined(mt.locus),
                                     (mt[cont_male_calls] / num_cont_males))
                               .when(mt.locus.in_x_nonpar() & hl.is_defined(mt.locus),
                                     (mt[cont_male_calls] + 2 * num_cont_females) /
                                     (num_cont_males + 2 * num_cont_females))
                               .default((mt[cont_male_calls] + mt[cont_female_calls]) /
                                        (num_cont_males + num_cont_females))),
               sexaware_cont_AC:
                              (hl.case()  # MINOR allele count
                               .when(mt.locus.in_y_nonpar(), mt[cont_male_homvars])
                               .when(mt.locus.in_x_nonpar(),
                                     mt[cont_male_homvars] + mt[cont_female_hets] + 2 * mt[cont_female_homvars])
                               .default(mt[cont_male_hets] + 2 * mt[cont_male_homvars] +
                                        mt[cont_female_hets] + 2 * mt[cont_female_homvars])),
               sexaware_cont_AN:
                              (hl.case()
                               .when(mt.locus.in_y_nonpar(), mt[cont_male_calls])
                               .when(mt.locus.in_x_nonpar(), mt[cont_male_calls] + 2 * mt[cont_female_calls])
                               .default(2 * mt[cont_male_calls] + 2 * mt[cont_female_calls]))})

        annotations_to_transfer.extend([cont_male_hets, cont_male_homvars, cont_male_calls, cont_female_hets,
                                        cont_female_homvars, cont_female_calls, sexaware_cont_call_rate,
                                        sexaware_cont_AC, sexaware_cont_AN])

    return mt, annotations_to_transfer


def annotate_variants_gnomad(mt, gnomad_ht):
    """
    Function to take variant information from Gnomad hail table and copy some (not all) annotations to matrix table
    rows.
    :param mt: matrix table to annotate
    :param gnomad_ht: string with file location + name of gnomad hail table to load
    :return: returns annotated matrix table
    """
    gnomad_sites = hl.read_table(gnomad_ht)

    mt = mt.annotate_rows(gnomad_freq=gnomad_sites.index(mt.row_key).freq,
                          gnomad_popmax=gnomad_sites.index(mt.row_key).popmax,
                          gnomad_faf=gnomad_sites.index(mt.row_key).faf,
                          gnomad_filters=gnomad_sites.index(mt.row_key).filters,
                          gnomad_rsid=gnomad_sites.index(mt.row_key).rsid,
                          gnomad_original_locus_37=gnomad_sites.index(mt.row_key).original_locus,
                          gnomad_original_alleles_37=gnomad_sites.index(mt.row_key).original_alleles)

    mt = mt.annotate_globals(gnomad_file=gnomad_ht)

    mt = mt.annotate_globals(gnomad_popmax_index_dict=gnomad_sites.index_globals().popmax_index_dict,
                             gnomad_freq_index_dict=gnomad_sites.index_globals().freq_index_dict,
                             gnomad_faf_index_dict=gnomad_sites.index_globals().faf_index_dict)

    return mt


def annotate_variants_gnomad_mismatch(mt, gnomad_mismatch_ht):
    """
    Imports a list of 'bad' variants that have significantly different frequencies between gnomad exomes and genomes
    in NFE population.
    :param mt: matrix table to annotate
    :param gnomad_mismatch_ht: string with file location + name of gnomad mismatch variant file to load
    :return: returns annotated matrix table
    """
    gnomad_mismatch_list = hl.read_table(gnomad_mismatch_ht)
    gnomad_mismatch_list = gnomad_mismatch_list.annotate(gnomad_mismatch=True)

    # gnomad_mismatch True/False boolean annotation
    mt = mt.annotate_rows(gnomad_mismatch_pvalue=gnomad_mismatch_list.index(mt.row_key).p_value,
                          gnomad_mismatch_variantid=gnomad_mismatch_list.index(mt.row_key).variant,
                          gnomad_mismatch=gnomad_mismatch_list.index(mt.row_key).gnomad_mismatch)

    # Fill in empty values for gnomad mismatch with False
    mt = mt.annotate_rows(gnomad_mismatch=hl.or_else(mt.gnomad_mismatch, False))

    mt = mt.annotate_globals(gnomad_mismatch_file=gnomad_mismatch_ht)

    return mt


def annotate_variants_cadd(mt, cadd_ht):
    """
    Function to take the file name of a hail table containing CADD information per variant, and annotate a given matrix
    table's variants with CADD information.

    :param mt: matrix table to annotatate
    :param cadd_ht: string with file name + location of hail table to read, containing CADD values.
    :return: returns annotated matrix table
    """
    cadd_variants = hl.read_table(cadd_ht)

    # Annotate variants with all CADD annotations, renaming them
    mt = mt.annotate_rows(CADD_RawScore=cadd_variants.index(mt.row_key).RawScore,
                          CADD_PHRED=cadd_variants.index(mt.row_key).PHRED)
    mt = mt.annotate_globals(cadd_file=cadd_ht)

    return mt


def annotate_variants_mpc(mt, mpc_ht):
    """
    Function to take MPC values and annotate rows of a matrix table with MPC values and information.

    :param mt: matrix table to annotate
    :param mpc_ht: string with file name + location of hail table to read, containing MPC values.
    :return: returns annotated matrix table
    """
    mpc_variants = hl.read_table(mpc_ht)

    # Annotate variants with all MPC annotations
    mt = mt.annotate_rows(**mpc_variants.index(mt.row_key))
    mt = mt.annotate_globals(mpc_file=mpc_ht)

    return mt


def test_missing_annotations(mt):
    """
    Function to test whether main annotations are missing and return tuple of boolean values for each.
    """
    try:
        test = mt.gnomad_freq
        gnomad_missing = False
    except:
        gnomad_missing = True

    try:
        test = mt.gnomad_mismatch_pvalue
        mismatch_missing = False
    except:
        mismatch_missing = True

    try:
        test = mt.CADD_PHRED
        cadd_missing = False
    except:
        cadd_missing = True

    try:
        test = mt.MPC
        mpc_missing = False
    except:
        mpc_missing = True

    return gnomad_missing, mismatch_missing, cadd_missing, mpc_missing
