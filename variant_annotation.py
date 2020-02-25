"""
Functions to annotate variants in a matrix table.

Created by Lea Urpa in 2019, some from Mitja's original code.
"""

import hail as hl


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
