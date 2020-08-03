"""
Functions to annotate variants in a matrix table.

Created by Lea Urpa in 2019, some from Mitja's original code.
"""

import hail as hl
import pandas as pd
import logging


def annotate_variants(mt):
    '''
    Takes matrix table and annotates variants with gene, LOF and missense annotations by parsing VEP annotations.

    :param mt: matrix table to annotate
    :return: returns matrix table with new row annotations gene, LOF, and missense.
    '''
    try:
        mt.row.was_split.describe()
    except Exception as e:
        print('Split multi-allelics before running!')
        print(e)
        return

    # If there is no canonical and protein-coding transcript consequence for that variant,
    # give the gene corresponding to the most severe consequence.
    # If there is a canonical and protein-coding transcript consequence for that variant,
    # give the gene symbol associated with that transcript consequence.
    canon_pc = mt.row.vep.transcript_consequences.filter(lambda x:
                                                         (x.canonical == 1) & (x.biotype == 'protein_coding'))
    most_severe = mt.row.vep.transcript_consequences.filter(lambda x:
                                                            x.consequence_terms.contains(
                                                                mt.row.vep.most_severe_consequence))

    mt = mt.annotate_rows(gene=hl.cond(hl.any(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding'),
                                       mt.row.vep.transcript_consequences),
                                       canon_pc.map(lambda x: x.gene_symbol),
                                       most_severe.map(lambda x: x.gene_symbol)))

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

    mt = mt.annotate_rows(LOF=hl.cond(canon_bool | non_canon_bool, True, False))

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

    mt = mt.annotate_rows(missense=hl.cond((canon_missense_bool | noncanon_missense_bool |
                                            canon_inframe_bool | noncanon_inframe_bool), True, False))

    # If the most severe consequence is "synonymous_variant", true else false
    mt = mt.annotate_rows(synonymous=hl.cond(mt.row.vep.most_severe_consequence == "synonymous_variant", True, False))

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


def annotate_cadd(mt, input_dir='.', cadd_filestem='dbNSFP4.0b2a_hg19CADD.chr'):
    '''
    Takes a matrix table keyed by locus/alleles (GrCH37/hg19) and annotates it with CADD information.

    :param mt: matrix table keyed by locus/alleles
    :param cadd_filestem: string referring to filestem for chrom-separated dbnsfp CADD info
    :param input_dir: string referring to input directory for dbnsfp data.
    :return: returns annotated matrix table with new row fields 'CADD_raw', 'CADD_raw_rankscore', and 'CADD_phred'
    '''

    # Run first for X chrom
    dbnsfp = hl.read_table(input_dir + cadd_filestem + 'X.ht')
    print('Annotating CADD for chr X')
    mt = mt.annotate_rows(CADD_raw=dbnsfp[mt.locus, mt.alleles].CADD_raw,
                          CADD_raw_rankscore=dbnsfp[mt.locus, mt.alleles].CADD_raw_rankscore,
                          CADD_phred=dbnsfp[mt.locus, mt.alleles].CADD_phred)

    # Run for the rest of the contigs
    iterators = [str(i).zfill(1) for i in range(1, 23)]
    iterators.extend(['Y', 'M'])

    for chrom in iterators:
        dbnsfp = hl.read_table(input_dir + cadd_filestem + str(chrom) + '.ht')
        print('Annotating CADD for chr ' + str(chrom))
        # If annotation is missing, input dbnsfp value
        mt = mt.annotate_rows(CADD_raw=hl.or_else(mt.CADD_raw, dbnsfp[mt.locus, mt.alleles].CADD_raw),
                              CADD_raw_rankscore=hl.or_else(mt.CADD_raw_rankscore,
                                                            dbnsfp[mt.locus, mt.alleles].CADD_raw_rankscore),
                              CADD_phred=hl.or_else(mt.CADD_phred, dbnsfp[mt.locus, mt.alleles].CADD_phred))

    mt = mt.annotate_globals(CADD_filestem=cadd_filestem)

    return mt


def annotate_mpc(mt, mpc):
    """
    Annotates an mt (keyed by locus/alleles) with missense badness- POLYPHEN2 - constraint (MPC) score.

    :param mt: matrix table, keyed by locus and alleles.
    :param mpc: hail table of mpc values for loci in the genome
    :return: returns annotated matrix table.
    """

    mt = mt.annotate_rows(MPC=mpc[mt.locus, mt.alleles].MPC,
                          mpc_ENST=mpc[mt.locus, mt.alleles].ENST,
                          mpc_ENSG=mpc[mt.locus, mt.alleles].ENSG,
                          mpc_gene_name=mpc[mt.locus, mt.alleles].gene_name,
                          mpc_consequence=mpc[mt.locus, mt.alleles].Consequence,
                          mpc_SIFT=mpc[mt.locus, mt.alleles].SIFT,
                          mpc_context=mpc[mt.locus, mt.alleles].context,
                          mpc_PolyPhen=mpc[mt.locus, mt.alleles].PolyPhen,
                          mpc_mis_badness=mpc[mt.locus, mt.alleles].mis_badness,
                          mpc_obs_exp=mpc[mt.locus, mt.alleles].obs_exp)

    return mt


def sex_aware_variant_annotations(mt, sex_col='sex', male_tag='Male', female_tag='Female'):
    '''
    Creates sex-aware variant annotations for call rate, allele count, and allele number.

    :param mt: matrix table to annotate
    :param sex_col: string referring to column in the matrix table giving sex information
    :param male_tag: string or boolean in the column referring to males
    :param female_tag: string or boolean in the column referring to females
    :return: Returns matrix table with new row annotations male_hets, male_homvars, male_calls, female_hets,
    female_homvars, female_calls, sexaware_call_rate, sexaware_ac and sexaware_an.
    '''
    num_males = mt.aggregate_cols(hl.agg.count_where(mt[sex_col] == male_tag))
    num_females = mt.aggregate_cols(hl.agg.count_where(mt[sex_col] == female_tag))

    mt = mt.annotate_rows(male_hets=hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) & (mt[sex_col] == male_tag)),
                          male_homvars=hl.agg.count_where(mt.GT.is_hom_var() & hl.is_defined(mt.GT) & (mt[sex_col] == male_tag)),
                          male_calls=hl.agg.count_where(hl.is_defined(mt.GT) & (mt[sex_col] == male_tag)),
                          female_hets=hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT) & (mt[sex_col] == female_tag)),
                          female_homvars=hl.agg.count_where(mt.GT.is_hom_var() & hl.is_defined(mt.GT) & (mt[sex_col] == female_tag)),
                          female_calls=hl.agg.count_where(hl.is_defined(mt.GT) & (mt[sex_col] == female_tag)))

    mt = mt.annotate_rows(sexaware_call_rate=(hl.case()
                                              .when(mt.locus.in_y_nonpar(), (mt.male_calls / num_males))
                                              .when(mt.locus.in_x_nonpar(),
                                                    (mt.male_calls + 2*mt.female_calls) / (num_males + 2*num_females))
                                              .default((mt.male_calls + mt.female_calls) / (num_males + num_females))),
                          sexaware_AC=(hl.case()  # MINOR allele count
                                       .when(mt.locus.in_y_nonpar(), mt.male_homvars)
                                       .when(mt.locus.in_x_nonpar(), mt.male_homvars + mt.female_hets + 2*mt.female_homvars)
                                       .default(mt.male_hets + 2*mt.male_homvars + mt.female_hets + 2*mt.female_homvars)),
                          sexaware_AN=(hl.case()
                                       .when(mt.locus.in_y_nonpar(), mt.male_calls)
                                       .when(mt.locus.in_x_nonpar(), mt.male_calls + 2*mt.female_calls)
                                       .default(2*mt.male_calls + 2*mt.female_calls)))
    logging.info('Completed sex-aware variant annotations.')

    return mt


def annotate_vars_gnomad_mismatch(mt, gnomad_filters_file):
    '''
    This function annotates the matrix table based on Elizabeth Atkins' list of variants
    with significantly different allele frequencies in Gnomad exomes and genomes
    within a population (NFE) which indicates that these variants are not sequenced well.

    :param mt: the matrix table to annotate
    :param gnomad_filters_file: a string containing the (google bucket, gs://) location of the
    tsv containing the chrom:pos:ref:alt variant, pvalue, odds ratio, genomes and exomes AFs.
    :return: Returns the annotated matrix table, with new row annotations gnomad_mismatch_pvalue,
    gnomad_mismatch_odds_ratio, gnomad_mismatch_genomes_AF_NFE and gnomad_mismatch_exomes_AF_NFE.
    '''

    # Load data
    gnomad_mismatch = hl.import_table(gnomad_filters_file, impute=True)

    # Convert to pandas to parse chrom:pos:ref:alt
    gnomad_mm_pd = gnomad_mismatch.to_pandas()
    loci = gnomad_mm_pd.V.str.split(pat=":", expand=True)
    loci['alleles'] = loci[[2, 3]].values.tolist()
    gnomad_mm_pd_loci = pd.concat([gnomad_mm_pd, loci], axis=1)

    # Convert back to hail table and create locus column, key table by these
    gnomad_mismatch_loci = hl.Table.from_pandas(gnomad_mm_pd_loci)
    gnomad_mismatch_loci = gnomad_mismatch_loci.annotate(
                            locus=hl.locus(gnomad_mismatch_loci['0'],
                                    hl.int(gnomad_mismatch_loci['1'])))
    gnomad_mismatch_loci = gnomad_mismatch_loci.key_by('locus','alleles')

    # Annotate
    mt = mt.annotate_rows(gnomad_mismatch_pvalue=gnomad_mismatch_loci[mt.locus, mt.alleles].p_value,
                          gnomad_mismatch_odds_ratio=gnomad_mismatch_loci[mt.locus, mt.alleles].odds_ratio,
                          gnomad_mismatch_genomes_AF_NFE=gnomad_mismatch_loci[mt.locus, mt.alleles].genomes_AF_NFE,
                          gnomad_mismatch_exomes_AF_NFE=gnomad_mismatch_loci[mt.locus, mt.alleles].exomes_AF_NFE)

    mt = mt.annotate_globals(gnomad_mismatch_file=gnomad_filters_file)

    return mt


def annotate_variants_gnomad(mt, gnomad_ht_file):
    '''
    Annotates variants based on gnomad information, including random forest filter values
    and population allele frequencies.

    :param mt: matrix table to annotate
    :param gnomad_ht_file: a string containing the (google bucket, gs://) location of the
    hail table to read.
    :return: returns the annotated matrix table, with new row annotations for gnomad_rf values,
    gnomad filter allele frequency (faf), rsid, population allele frequencies, and popmax.
    '''
    gnomad_sites = hl.read_table(gnomad_ht_file)

    # Annotate our data with Gnomad info
    mt = mt.annotate_rows(gnomad_filters=gnomad_sites[mt.locus, mt.alleles].filters,  # this is the one we want
                          gnomad_faf=gnomad_sites[mt.locus, mt.alleles].faf,
                          gnomad_rsid=gnomad_sites[mt.locus, mt.alleles].rsid,
                          # gnomad_afs=gnomad_sites[mt.locus, mt.alleles].vep.colocated_variants,
                          gnomad_popmax=gnomad_sites[mt.locus, mt.alleles].popmax,
                          gnomad_freqs=gnomad_sites[mt.locus, mt.alleles].freq)

    # Hard filtering was for Gnomad's QC measures (what I call analytical variant qc)
    mt = mt.annotate_rows(gnomad_fail_hard_filters=gnomad_sites[mt.locus, mt.alleles].fail_hard_filters)

    mt = mt.annotate_globals(gnomad_file=gnomad_ht_file,
                             gnomad_freq_index_dict=gnomad_sites.index_globals().freq_index_dict,
                             gnomad_popmax_index_dict=gnomad_sites.index_globals().popmax_index_dict,
                             gnomad_faf_index_dict=gnomad_sites.index_globals().faf_index_dict)

    return mt

