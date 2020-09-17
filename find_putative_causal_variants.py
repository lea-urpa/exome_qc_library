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

    mt = mt.checkpoint(args.output_stem + "_non_monomorphic_tmp.mt")

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
    ############################################################
    # Get count of samples that are cases and controls, report #
    ############################################################
    logging.info("Annotating het and hom var carrier counts in controls to variants.")

    case_count = mt.aggregate_cols(hl.agg.count_where(mt[args.pheno_col] == True))
    control_count = mt.aggregate_cols(hl.agg.count_where(mt[args.pheno_col] == False))
    missing = mt.aggregate_cols(hl.agg.count_where(~hl.is_defined(mt[args.pheno_col])))

    logging.info(f"Number of controls in dataset: {control_count}")
    logging.info(f"Number of cases in dataset: {case_count}")
    logging.info(f"Samples missing case/control information: {missing}")
    if missing > 0:
        logging.info(f"Warning- samples missing case/control status will be generally ignored in this pipeline.")

    ##################################################
    # Annotate control/case het count + homvar count #
    ##################################################
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

    ################################################################################
    # Count control homvar females and het males for hemizygous mutations on X chr #
    ################################################################################
    mt = mt.annotate_rows(control_homvar_count_female=
                          hl.agg.filter((mt[args.pheno_col] == False) & hl.is_defined(mt[args.pheno_col]) &
                                        (mt[args.female_col] == True) & hl.is_defined(mt[args.female_col]),
                                        hl.agg.count_where(mt.GT.is_hom_var())))

    mt = mt.annotate_rows(control_het_count_male=
                          hl.agg.filter((mt[args.pheno_col] == False) & hl.is_defined(mt[args.pheno_col]) &
                                        (mt[args.female_col] == False) & hl.is_defined(mt[args.female_col]),
                                        hl.agg.count_where(mt.GT.is_het())))

    mt = mt.annotate_rows(case_homvar_count_female=
                          hl.agg.filter((mt[args.pheno_col] == True) & hl.is_defined(mt[args.pheno_col]) &
                                        (mt[args.female_col] == True) & hl.is_defined(mt[args.female_col]),
                                        hl.agg.count_where(mt.GT.is_hom_var())))

    mt = mt.annotate_rows(case_het_count_male=
                          hl.agg.filter((mt[args.pheno_col] == True) & hl.is_defined(mt[args.pheno_col]) &
                                        (mt[args.female_col] == False) & hl.is_defined(mt[args.female_col]),
                                        hl.agg.count_where(mt.GT.is_het())))

    return mt


def annotate_variants(mt, args):
    """
    Annotates matrix table variants with CADD, MPC and gnomad info.
    :param mt: matrix table to annotate
    :param args: arguments for cadd, mpc, and gnomad hail table locations
    :return: returns annotated matrix table
    """
    logging.info("Annotating matrix table with CADD, MPC and Gnomad.")
    checkpoint_name = args.output_stem + "_annotation_tmp.mt"

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

    mt = mt.checkpoint(checkpoint_name, overwrite=True)

    return mt


def annotate_population_thresholds(mt, args):
    """
    Annotate with True/False categories based on population parameters for recessive and dominant inheritance models
    in Gnomad population counts and allele frequencies and in control counts in dataset.

    :param mt: matrix table to annotate
    :param args: arguments giving max allowed carrier counts, max allele frequencies, and gnomad population index
    :return: returns annotated matrix table
    """
    logging.info("Annotating with boolean columns for whether variant fulfills population threshold criteria.")
    checkpoint_name = args.output_stem + "_annotation_tmp.mt"

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
                          dominant_rare_controls=rows[mt.locus, mt.alleles].hemizygous_rare_gnomad,
                          recessive_rare_controls=rows[mt.locus, mt.alleles].recessive_rare_controls,
                          hemizygous_rare_controls=rows[mt.locus, mt.alleles].hemizygous_rare_controls)

    mt = mt.checkpoint(checkpoint_name, overwrite=True)

    return mt


def annotate_genes(mt, args):
    """
    Annotates gene set information, if disease geneset given.

    :param mt: Matrix table to annotate gene set information to
    :param args: arguments giving gene list, disease gene sep, gene col name, and output stem
    :return: Returns matrix table with rows annotated with T/F column of whether variant in gene set of interest,
    and column of allelic requirement strings.
    """
    ###########################################################################
    # Pull MT rows, select only locus, alleles, and gene, checkpoint, explode #
    ###########################################################################
    rows = mt.rows()
    rows = rows.key_by()
    rows = rows.select("locus", "alleles", "gene")
    rows = rows.checkpoint(args.output_stem + "_rows_tmp.ht")
    rows = rows.explode(rows.gene)
    rows = rows.key_by(rows.gene)

    if args.gene_list is not None:
        logging.info("Annotating dataset with disease gene information.")
        gene_table = hl.import_table(args.gene_list, delimiter=args.disease_gene_sep)
        gene_table = gene_table.transmute(gene=gene_table[args.gene_col_name])
        gene_table = gene_table.key_by('gene')

        ##################################################################################
        # Add disease gene information as column, split allelic requirements and explode #
        ##################################################################################
        rows = rows.annotate(**gene_table[rows.gene])
        rows = rows.annotate(**{args.gene_set_name: hl.cond(hl.is_defined(rows[args.allelic_requirement_col]), True, False)})

        rows = rows.annotate(allelic_requirement=rows[args.allelic_requirement_col].split(","))
        rows = rows.explode(rows.allelic_requirement)

        ######################################################
        # Annotate inheritance based on allelic requirements #
        ######################################################
        dominant_terms = ['biallelic', 'uncertain', 'digenic']
        recessive_terms = ['monoallelic', 'imprinted', 'x-linked dominant', 'x-linked over-dominance', 'uncertain',
                           'digenic', 'mosaic']
        rows = rows.annotate(inheritance=hl.case()
                             .when(hl.array(dominant_terms).contains(rows.allelic_requirement), 'dominant')
                             .when(hl.array(recessive_terms).contains(rows.allelic_requirement), 'recessive')
                             .when(rows.allelic_requirement == 'hemizygous', 'hemizygous')
                             .or_missing())
        ###################################
        # Group table by gene, checkpoint #
        ###################################
        disease_genes = rows.group_by('locus', 'alleles').aggregate(
                            **{'gene': hl.array(hl.agg.collect_as_set(rows.gene)),
                               args.gene_set_name: hl.array(hl.agg.collect_as_set(rows[args.gene_set_name])),
                               'allelic_requirement': hl.array(hl.agg.collect_as_set(rows.allelic_requirement)),
                               'inheritance': hl.array(hl.agg.collect_as_set(rows.inheritance))})
        disease_genes = disease_genes.key_by(disease_genes.locus, disease_genes.alleles)

        #######################################
        # Annotate rows with disease genes ht #
        #######################################
        mt = mt.annotate_rows(allelic_requirement=disease_genes[mt.locus, mt.alleles].allelic_requirement)
        mt = mt.annotate_rows(**{args.gene_set_name: disease_genes[mt.locus, mt.alleles][args.gene_set_name]})
        mt = mt.annotate_rows(inheritance=disease_genes[mt.locus, mt.alleles].inheritance)
    else:
        mt = mt.annotate_rows(allelic_requirement=hl.empty_array(hl.tstr))
        mt = mt.annotate_rows(**{args.gene_set_name: hl.empty_array(hl.tset)})
        mt = mt.annotate_rows(inheritance=hl.empty_array(hl.tstr))

    if args.gnomad_gene_metrics is not None:
        logging.info("Annotating genes with pLI metrics.")
        gene_metrics = hl.import_table(args.gnomad_gene_metrics, types={'pLI': hl.tfloat64}, key='gene')

        #################################
        # Add pLI information as column #
        #################################
        rows = rows.annotate(pLI=gene_metrics[rows.gene].pLI)

        ###################################
        # Group table by gene, checkpoint #
        ###################################
        gene_metrics = rows.group_by('locus', 'alleles').aggregate(pLI=hl.agg.collect(rows.pLI))
        gene_metrics = gene_metrics.key_by(gene_metrics.locus, gene_metrics.alleles)

        ###########################################
        # Annotate matrix table with gene metrics #
        ###########################################
        mt = mt.annotate_rows(pLI=gene_metrics[mt.locus, mt.alleles].pLI)
    else:
        mt = mt.annotate_rows(pLI=hl.empty_array(hl.tstr))

    return mt


def find_putative_causal_variants(mt, args):
    """
    Using rules for inheritance type, gene set, and given mutation consequences, searches for variants that satisfy
    criteria for gene set inclusion, predicted functional consequences, presence in controls + gnomad

    :param mt: matrix table to look for putative causal variants
    :param args: arguments for output info and thresholds for presence in controls and missense damaging scores
    :return: returns table of putative causal variants with carriers annotated
    """
    ############################
    # Pull rows and checkpoint #
    ############################
    rows = mt.rows()
    rows = rows.checkpoint(args.output_stem+ "_rows_tmp.ht", overwrite=True)

    #################################
    # Annotate putative causal: LOF #
    #################################
    rows = rows.annotate(putative_causal=hl.empty_array(hl.tstruct))

    inheritances = ['dominant', 'recessive', 'hemizygous']
    gene_sets = [args.gene_set_name, 'high_pLI']
    consequences = ['LOF', 'missense', 'noncoding']

    for consequence in consequences:
        for gene_set in gene_sets:
            for inheritance in inheritances:
                annotation = {}  # reset annotation for each combination of consequence, gene set, inheritance

                # These must all be done in innermost loop, otherwise you reference 2 different HT objects

                ###############################################
                # Create boolean for the mutation consequence #
                ###############################################
                if consequence is 'LOF':
                    consequence_bool = (rows.LOF == True)  # variant is LOF
                    annotation['consequence'] = 'lof'
                elif consequence is 'missense':
                    annotation['consequence'] = 'missense'
                    if inheritance is 'dominant':
                        # missense, high MPC score
                        consequence_bool = (rows.missense == True) & (rows.MPC > args.mpc_cutoff)
                    if inheritance is 'recessive':
                        # missense, high CADD score
                        consequence_bool = (rows.missense == True) & (rows.CADD_PHRED > args.cadd_cutoff)
                    if inheritance is 'hemizygous':
                        # missense, either high CADD score or high MPC score
                        consequence_bool = ((rows.missense == True) &
                                            ((rows.CADD_PHRED > args.cadd_cutoff) | (rows.MPC > args.mpc_cutoff)))
                elif consequence is 'noncoding':  # noncoding
                    consequence_bool = (rows.LOF == False) & (rows.missense == False)  # neither LOF nor missense
                    annotation['consequence'] = 'noncoding'

                ################################
                # Create boolean for gene sets #
                ################################
                if gene_set is args.gene_set_name:
                    gene_set_bool = hl.any(lambda x: x == True, rows[args.gene_set_name])  # variant is in gene set
                    annotation['gene_set'] = args.gene_set_name
                elif gene_set is 'high_pLI':
                    gene_set_bool = (
                            hl.any(lambda x: x > 0.9, rows.pLI) & hl.is_defined(rows.pLI) &  # variant is high pLI
                            ~hl.any(lambda x: x == True, rows[args.gene_set_name]))  # variant is not in gene set
                    annotation['gene_set'] = 'high_pLI'

                ###################################
                # Create boolean for inheritances #
                ###################################
                if inheritance is 'dominant':
                    annotation['inheritance'] = 'dominant'
                    if gene_set is args.gene_set_name:
                        inheritance_bool = (hl.any(lambda x: x == 'dominant', rows.inheritance) &  # dom. evidence
                                            (rows.dominant_rare_controls == True) &  # rare in controls, dominant
                                            ((rows.dominant_rare_gnomad == True) |
                                             ~hl.is_defined(rows.dominant_rare_gnomad)) &  # rare in gnomad
                                            (rows.case_het_count > 0))  # at least one het case
                    else:
                        inheritance_bool = ((rows.dominant_rare_controls == True) &  # rare in controls, dominant
                                            ((rows.dominant_rare_gnomad == True) |
                                             ~hl.is_defined(rows.dominant_rare_gnomad)) &  # rare in gnomad
                                            (rows.case_het_count > 0))  # at least one het case
                if inheritance is 'recessive':
                    annotation['inheritance'] = 'recessive'
                    if gene_set is args.gene_set_name:
                        inheritance_bool = (hl.any(lambda x: x == 'recessive', rows.inheritance) &  # rec. evidence
                                            (rows.recessive_rare_controls == True) &  # rare in controls, recessive
                                            ((rows.recessive_rare_gnomad == True) |
                                             ~hl.is_defined(rows.recessive_rare_gnomad)) &  # rare in gnomad
                                            (rows.case_homvar_count > 0))  # at least one homvar case
                    else:
                        inheritance_bool = ((rows.locus.contig != 'chrX') &  # locus not on chromosome x
                                            (rows.recessive_rare_controls == True) &  # rare in controls, recessive
                                            ((rows.recessive_rare_gnomad == True) |
                                             ~hl.is_defined(rows.recessive_rare_gnomad)) &  # rare in gnomad
                                            (rows.case_homvar_count > 0))  # at least one homvar case

                if inheritance is 'hemizygous':
                    annotation['inheritance'] = 'hemizygous'
                    if gene_set is args.gene_set_name:
                        inheritance_bool = (hl.any(lambda x: x == 'hemizygous', rows.inheritance) &  # hemi. evidence
                                            (rows.hemizygous_rare_controls == True) &  # rare in controls, recessive
                                            ((rows.hemizygous_rare_gnomad == True) |
                                             ~hl.is_defined(rows.hemizygous_rare_gnomad)) &  # rare in gnomad
                                            ((rows.case_homvar_count_female > 0) |
                                             (rows.case_het_count_male > 0)))  # one male het or female homvar
                    else:
                        inheritance_bool = ((rows.hemizygous_rare_controls == True) &  # rare in controls
                                            ((rows.hemizygous_rare_gnomad == True) |
                                             ~hl.is_defined(rows.hemizygous_rare_gnomad)) &  # rare in gnomad
                                            ((rows.case_homvar_count_female > 0) |
                                             (rows.case_het_count_male > 0)))  # one male het or female homvar

                ###########################################################################
                # Query each variant for consequence, gene set, and inheritance > causal? #
                ###########################################################################
                rows = rows.annotate(putative_causal=hl.cond(
                    consequence_bool & gene_set_bool & inheritance_bool,
                    rows.putative_causal.append(annotation), rows.putative_causal))

    ##################################################################
    # Annotate putative causal to matrix table, filter, get carriers #
    ##################################################################
    mt = mt.annotate_rows(putative_causal=rows[mt.locus, mt.alleles].putative_causal)
    mt = mt.filter_rows(hl.len(mt.putative_causal) > 0, keep=True)

    # Annotate carriers now that there should be not too many
    mt = mt.annotate_rows(het_carriers=hl.agg.filter(mt.GT.is_het() & hl.is_defined(mt.GT), hl.agg.collect(mt.s)),
                          homvar_carriers=hl.agg.filter(mt.GT.is_hom_var() & hl.is_defined(mt.GT), hl.agg.collect(mt.s)),
                          hemizygous_carriers=hl.cond(
                              mt[args.female_col] == True, # If female, collect homvar carriers, else collect hets
                              hl.agg.filter(mt.GT.is_hom_var() & hl.is_defined(mt.GT), hl.agg.collect(mt.s)),
                              hl.agg.filter(mt.GT.is_het() & hl.is_defined(mt.GT), hl.agg.collect(mt.s))))

    mt = mt.write(args.output_stem + "_putative_causal_final.mt", overwrite=True)

    rows = mt.rows()
    rows = rows.checkpoint(args.output_stem + "_rows_tmp.ht", overwrite=True)
    counter = rows.aggregate(hl.agg.counter(hl.len(rows.putative_causal)))
    logging.info(f"Counter of number of putative causal annotations per variant: {counter}")

    return rows


def annotate_denovos(rows, args):

    #############################################################################################
    # Pull dominant, recessive and hemizgous to separate tables, explode by appropriate carrier #
    #############################################################################################
    dominant = rows.filter(hl.any(lambda x: x['inheritance'] == 'dominant', rows.putative_causal), keep=True)
    dominant = dominant.drop('homvar_carriers', 'hemizygous_carriers')
    dominant = dominant.explode(dominant.het_carriers)
    dominant = dominant.transmute(id=dominant.het_carriers)

    recessive = rows.filter(hl.any(lambda x: x['inheritance'] == 'recessive', rows.putative_causal), keep=True)
    recessive = recessive.drop('het_carriers', 'hemizygous_carriers')
    recessive = recessive.explode(recessive.homvar_carriers)
    recessive = recessive.transmute(id=recessive.homvar_carriers)

    hemizygous = rows.filter(hl.any(lambda x: x['inheritance'] == 'hemizygous', rows.putative_causal), keep=True)
    hemizygous = hemizygous.drop('het_carriers', 'homvar_carriers')
    hemizygous = hemizygous.explode('hemizygous_carriers')
    hemizygous = hemizygous.transmute(id=hemizygous.hemizygous_carriers)

    ##############################
    # Annotate de novo mutations #
    ##############################
    if args.denovo_ht is not None:
        de_novos = hl.read_table(args.denovo_ht)

        dominant = dominant.key_by('locus', 'alleles', 'id')
        dominant = dominant.annotate(p_de_novo=de_novos[de_novos.locus, de_novos.alleles, de_novos.id].p_de_novo,
                                     denovo_confidence=de_novos[de_novos.locus, de_novos.alleles, de_novos.id].confidence)

        recessive = recessive.key_by('locus', 'alleles', 'id')
        recessive = recessive.annotate(p_de_novo=de_novos[de_novos.locus, de_novos.alleles, de_novos.id].p_de_novo,
                                       denovo_confidence=de_novos[de_novos.locus, de_novos.alleles, de_novos.id].confidence)

        hemizygous = hemizygous.key_by('locus', 'alleles', 'id')
        hemizygous = hemizygous.annotate(p_de_novo=de_novos[de_novos.locus, de_novos.alleles, de_novos.id].p_de_novo,
                                         denovo_confidence=de_novos[de_novos.locus, de_novos.alleles, de_novos.id].confidence)

    causal_vars = dominant.union(recessive, hemizygous, unify=True)

    if args.denovo_ht is None:
        causal_vars = causal_vars.annotate(p_de_novo=hl.null(hl.tfloat64), denovo_confidence=hl.null(hl.tstr))

    ##############################
    # Annotate evidence category #
    ##############################
    causal_vars = causal_vars.annotate(de_novo=hl.cond(
        hl.array(['HIGH', 'MEDIUM']).contains(causal_vars.denovo_confidence) &
        hl.is_defined(causal_vars.denovo_confidence),
        True, False))

    # Start with weakly suggestive, stronger evidence overwrites if variant on multiple genes
    causal_vars = causal_vars.annotate(
        putative_category=hl.case()
        .when(hl.any(lambda x: (x['consequence'] == 'noncoding') & (x['gene_set'] == args.gene_set_name) &
                               (causal_vars.de_novo == True),
                     causal_vars.putative_causal), 'weakly_suggestive')
        .when(hl.any(lambda x: (x['consequence'] == 'missense') & (x['gene_set'] == args.gene_set_name) &
                               (causal_vars.de_novo == False),
                     causal_vars.putative_causal), 'weakly_suggestive')
        .when(hl.any(lambda x: (x['consequence'] == 'lof') & (x['gene_set'] == 'high_pLI') &
                               (causal_vars.de_novo == False),
                     causal_vars.putative_causal), 'weakly_suggestive')
        .when(hl.any(lambda x: (x['consequence'] == 'missense') & (x['gene_set'] == args.gene_set_name) &
                               (causal_vars.de_novo == True),
                     causal_vars.putative_causal), 'moderately_suggestive')
        .when(hl.any(lambda x: (x['consequence'] == 'lof') & (x['gene_set'] == 'high_pLI') &
                               (causal_vars.de_novo == True),
                     causal_vars.putative_causal), 'moderately_suggestive')
        .when(hl.any(lambda x: (x['consequence'] == 'lof') & (x['gene_set'] == args.gene_set_name),
                     causal_vars.putative_causal), 'strongly_suggestive')
        )

    counter = causal_vars.aggregate(hl.agg.counter(causal_vars.putative_category))
    logging.info(f"Number of causal variants per individual in each putative category: {counter}")

    ###########################
    # Export table + tsv file #
    ###########################
    causal_vars.write(args.output_dir + "_causal_vars_table_final.ht", overwrite=True)
    causal_vars.export(args.output_dir + "_causal_vars_table_final.txt")


if __name__ == "__main__":
    import hail as hl
    hl.init()

    ###################
    # Parse arguments #
    ###################
    parser = argparse.ArgumentParser(description="v9 exome sequencing dataset quality control pipeline.")
    parser.add_argument("-mt", type=str, help="Input matrix table to run analysis on")
    parser.add_argument("--pheno_col", required=True, type=str, help="Col annotation giving case status, true/false.")
    parser.add_argument("--female_col", type=str, help="Col annotation giving female status, true/false.",
                        default="is_female_imputed")
    parser.add_argument("--de_novo_ht", type=str, help="Hail table output from hl.de_novo")
    parser.add_argument("--disease_genes", type=str, help="Name of file containing genes implicated in the disease.")
    parser.add_argument("--gene_col_name", type=str, default="gene.symbol")
    parser.add_argument("--allelic_requirement_col", type=str, default="allelic.requirement.final")
    parser.add_argument("--disease_gene_sep", type=str, default="\t")
    parser.add_argument("--gene_set_name", type=str, help='Name of the disease gene set for annotation.')
    parser.add_argument("--mpc_cutoff", type=float, default=2,
                        help="Threshold for damaging missense variants by MPC score, for dominant model.")
    parser.add_argument("--cadd_cutoff", type=float, default=20,
                        help="Threshold for damaging missense variants by CADD score, for recessive model.")
    parser.add_argument("--max_allowed_carrier_dominant", type=int, default=5,
                        help="Maximum allowed number of carriers for a variant for dominant inheritance model.")
    parser.add_argument("--max_allowed_homozygotes_recessive", type=int, default=5, help="Maximum allowed number of "
                        "homozygous carriers for a variant in recessive inheritance model.")
    parser.add_argument("--gnomad_AF_cutoff_recessive", type=float, default=0.01,
                        help="Max allele frequency in Gnomad to consider a variant damaging in recessive model.")
    parser.add_argument("--cadd_ht", required=True, type=str, help="File name of CADD hail table.")
    parser.add_argument("--mpc_ht", required=True, type=str, help="File name of MPC hail table.")
    parser.add_argument("--gnomad_ht", required=True, type=str, help="File name of gnomad hail table.")
    parser.add_argument("--gnomad_population", type=str, default="controls",
                        choices=['gnomad', 'non_topmed', 'non_neuro', 'non_cancer', 'controls'])
    parser.add_argument("--gnomad_mismatch_ht", type=str,
                        help="Hail table with list of Gnomad variants with significantly different AF between genomes"
                             "and exomes.")
    parser.add_argument("--gnomad_gene_metrics", type=str, help='text file with gnomad gene metrics')
    parser.add_argument("--output_name", required=True, type=str, help="Output name for files.")
    parser.add_argument("--output_dir", required=True, type=str, help="Output directory for output files.")
    parser.add_argument("--scripts_dir", required=True, type=str, help="Directory with scripts from this library")
    parser.add_argument("--log_dir", required=True, type=str, help="Output directory for logs.")
    args = parser.parse_args()

    args.output_stem = os.path.join(args.output_dir, args.output_name)
    if args.gene_set_name is None:
        args.gene_set_name = 'no_geneset'

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

    ##################################################
    # Annotate matrix table with various annotations #
    ##################################################
    var_mt = annotate_control_carriers(var_mt, args)
    var_mt = annotate_variants(var_mt, args)
    var_mt = annotate_population_thresholds(var_mt, args)
    var_mt = annotate_genes(var_mt, args)

    ##########################################
    # Run analysis to find putative variants #
    ##########################################
    variants = find_putative_causal_variants(var_mt, args)
    annotate_denovos(variants, args)

    ###########################
    # Copy logs and shut down #
    ###########################
    logging.info('Pipeline ran successfully! Copying logs and shutting down cluster in 10 minutes.')
    h.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.log_dir)
