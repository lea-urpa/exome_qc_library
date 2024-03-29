"""
Script to use predicted variant consequence and population information to find putative causal variants for cases.
"""
import time
import os
import logging
import argparse
import sys
import causal_variant_functions as cv
import samples_qc as sq
import variant_qc as vq
import subprocess
import hail as hl


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
    rows = rows.checkpoint(args.output_stem + "_rows_tmp3.ht", overwrite=True)

    #################################
    # Annotate putative causal: LOF #
    #################################
    rows = rows.annotate(putative_causal=hl.empty_array(hl.tdict(hl.tstr, hl.tstr)))

    inheritances = ['dominant', 'recessive', 'hemizygous']
    gene_sets = [args.gene_set_name]
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
                if gene_set is 'high_pLI':  # If no gene set given, look for variants in high pLI genes
                    gene_set_bool = rows.high_pLI == True # variant is high pLI
                    annotation['gene_set'] = 'high_pLI'
                else:
                    gene_set_bool = rows[args.gene_set_name] == True  # variant is in gene set
                    annotation['gene_set'] = args.gene_set_name

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
                query = rows.aggregate(hl.agg.count_where(consequence_bool & gene_set_bool & inheritance_bool))
                logging.info(f"Number of variants that are {consequence} in {gene_set} genes with "
                             f"{inheritance} inheritance: {query}")

                rows = rows.annotate(putative_causal=hl.cond(
                    consequence_bool & gene_set_bool & inheritance_bool,
                    rows.putative_causal.append(annotation), rows.putative_causal))

    ##################################################################
    # Annotate putative causal to matrix table, filter, get carriers #
    ##################################################################
    mt = mt.annotate_rows(putative_causal=rows[mt.locus, mt.alleles].putative_causal,
                          dominant_rare_gnomad=rows[mt.locus, mt.alleles].dominant_rare_gnomad,
                          recessive_rare_gnomad=rows[mt.locus, mt.alleles].recessive_rare_gnomad,
                          hemizygous_rare_gnomad=rows[mt.locus, mt.alleles].hemizygous_rare_gnomad,
                          dominant_rare_controls=rows[mt.locus, mt.alleles].dominant_rare_controls,
                          recessive_rare_controls=rows[mt.locus, mt.alleles].recessive_rare_controls,
                          hemizygous_rare_controls=rows[mt.locus, mt.alleles].hemizygous_rare_controls)

    mt.write(args.output_stem + "_putative_causal_allvariants.mt", overwrite=True)

    mt = mt.filter_rows(hl.len(mt.putative_causal) > 0, keep=True)

    # Annotate carriers now that there should be not too many
    mt = mt.annotate_rows(het_carriers=hl.agg.filter(mt.GT.is_het() & hl.is_defined(mt.GT), hl.agg.collect(mt.s)),
                          homvar_carriers=hl.agg.filter(mt.GT.is_hom_var() & hl.is_defined(mt.GT), hl.agg.collect(mt.s)),
                          het_male_carriers=hl.agg.filter(
                              mt.GT.is_het() & hl.is_defined(mt.GT) & (mt[args.female_col] == False),
                              hl.agg.collect(mt.s)),
                          homvar_female_carriers=hl.agg.filter(
                              mt.GT.is_hom_var() & hl.is_defined(mt.GT) & (mt[args.female_col] == True),
                              hl.agg.collect(mt.s)))

    mt = mt.annotate_rows(hemizygous_carriers=hl.flatten([mt.het_male_carriers, mt.homvar_female_carriers]))
    mt.write(args.output_stem + "_putative_causal_final.mt", overwrite=True)

    rows = mt.rows()
    rows = rows.checkpoint(args.output_stem + "_rows_tmp4.ht", overwrite=True)
    counter = rows.aggregate(hl.agg.counter(hl.len(rows.putative_causal)))
    logging.info(f"Counter of number of putative causal annotations per variant: {counter}")

    return rows, mt


def annotate_denovos_genotypes(rows, mt, args):

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
    if args.de_novo_ht is not None:
        de_novos = hl.read_table(args.de_novo_ht)
        de_novos = de_novos.key_by('locus', 'alleles', 'id')

        dominant = dominant.key_by('locus', 'alleles', 'id')
        dominant = dominant.annotate(
            p_de_novo=de_novos[dominant.locus, dominant.alleles, dominant.id].p_de_novo,
            denovo_confidence=de_novos[dominant.locus, dominant.alleles, dominant.id].confidence)

        recessive = recessive.key_by('locus', 'alleles', 'id')
        recessive = recessive.annotate\
            (p_de_novo=de_novos[recessive.locus, recessive.alleles, recessive.id].p_de_novo,
             denovo_confidence=de_novos[recessive.locus, recessive.alleles, recessive.id].confidence)

        hemizygous = hemizygous.key_by('locus', 'alleles', 'id')
        hemizygous = hemizygous.annotate(
            p_de_novo=de_novos[hemizygous.locus, hemizygous.alleles, hemizygous.id].p_de_novo,
            denovo_confidence=de_novos[hemizygous.locus, hemizygous.alleles, hemizygous.id].confidence)

    causal_vars = dominant.union(recessive, hemizygous, unify=True)
    causal_vars = causal_vars.checkpoint(args.output_stem + "_causalvars_tmp1.ht", overwrite=True)
    h.add_preemptibles(args.cluster_name, args.num_preemptibles)

    if args.de_novo_ht is None:
        causal_vars = causal_vars.annotate(p_de_novo=hl.null(hl.tfloat64), denovo_confidence=hl.null(hl.tstr))

    #####################################################
    # Annotate with entry information from matrix table #
    #####################################################
    entries = mt.entries()
    causal_vars = causal_vars.key_by()
    causal_vars = causal_vars.transmute(s=causal_vars.id)
    causal_vars = causal_vars.key_by('locus', 'alleles', 's')
    causal_vars = causal_vars.checkpoint(args.output_stem + "_causalvars_tmp2.ht", overwrite=True)

    causal_vars = causal_vars.annotate(
    GT=entries[causal_vars.locus, causal_vars.alleles, causal_vars.s].GT,
    final_failing_depth_quality=entries[causal_vars.locus, causal_vars.alleles, causal_vars.s].final_failing_depth_quality,
    final_failing_ab=entries[causal_vars.locus, causal_vars.alleles, causal_vars.s].final_failing_ab,
    AD=entries[causal_vars.locus, causal_vars.alleles, causal_vars.s].AD,
    DP=entries[causal_vars.locus, causal_vars.alleles, causal_vars.s].DP,
    GQ=entries[causal_vars.locus, causal_vars.alleles, causal_vars.s].GQ)

    ##############################
    # Annotate evidence category #
    ##############################
    causal_vars = causal_vars.annotate(de_novo=hl.cond(
        hl.array(['HIGH', 'MEDIUM']).contains(causal_vars.denovo_confidence) &
        hl.is_defined(causal_vars.denovo_confidence),
        True, False))

    # Start with weakly suggestive, stronger evidence overwrites if variant on multiple genes
    causal_vars = causal_vars.annotate(
        putative_category=(hl.case()
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
        .or_missing()))

    causal_vars = causal_vars.checkpoint(args.output_stem + "_causalvars_tmp3.ht", overwrite=True)

    counter = causal_vars.aggregate(hl.agg.counter(causal_vars.putative_category))
    logging.info(f"Number of causal variants per individual in each putative category: {counter}")

    ###########################
    # Export table + tsv file #
    ###########################
    causal_vars.write(args.output_stem + "_causal_vars_table_final.ht", overwrite=True)

    # Pull out relevant gnomad population information, drop large structs, and flatten
    causal_vars = causal_vars.annotate(
        popmax_population=causal_vars.gnomad_popmax[causal_vars.gnomad_popmax_index_dict['gnomad']].pop)
    causal_vars = causal_vars.annotate(
        gnomad_fin_freq=causal_vars.gnomad_freq[causal_vars.gnomad_freq_index_dict['gnomad_fin']],
        gnomad_controls_fin_freq=causal_vars.gnomad_freq[causal_vars.gnomad_freq_index_dict['controls_fin']],
        gnomad_nonneuro_fin_freq=causal_vars.gnomad_freq[causal_vars.gnomad_freq_index_dict['non_neuro_fin']],
        gnomad_popmax_gnomad=causal_vars.gnomad_popmax[causal_vars.gnomad_popmax_index_dict['gnomad']],
        gnomad_popmax_controls=causal_vars.gnomad_popmax[causal_vars.gnomad_popmax_index_dict['controls']],
        gnomad_popmax_nonneuro=causal_vars.gnomad_popmax[causal_vars.gnomad_popmax_index_dict['non_neuro']])
    causal_vars = causal_vars.drop('gnomad_freq', 'gnomad_popmax', 'vep')
    causal_vars = causal_vars.flatten()

    causal_vars.export(args.output_stem + "_causal_vars_table_final.txt")


if __name__ == "__main__":
    import hail as hl
    hl.init()

    ###################
    # Parse arguments #
    ###################
    parser = argparse.ArgumentParser(description="Exome sequencing dataset quality control pipeline.")
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
    parser.add_argument("--pLI_cutoff", type=float, help="Cutoff for high pLI genes", default=0.9)
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
    parser.add_argument("--cluster_name", required=True, type=str, help="Cluster name for adding preemptibles.")
    parser.add_argument("--num_preemptibles", required=True, type=int, help="Number of preemptible nodes to add.")
    args = parser.parse_args()

    args.output_stem = os.path.join(args.output_dir, args.output_name)
    if args.disease_genes is None:
        args.gene_set_name = 'high_pLI'

    if (args.disease_genes is None) and (args.gnomad_gene_metrics is None):
        logging.error("Error! If not giving disease gene list, gnomad gene metrics to get pLI values must be given.")
        exit()

    if (args.disease_genes is not None) and (args.gene_set_name is None):
        logging.error("Error! If giving disease gene list, --gene_set_name must also be given.")
        exit()

    scripts = ["variant_annotation.py", "helper_scripts.py", "utils.py"]
    for script in scripts:
        hl.spark_context().addPyFile(os.path.join(args.scripts_dir, script))

    import variant_annotation as va
    import helper_scripts as h
    import utils

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
    utils.add_secondary(args.cluster_name, args.num_secondary, region=args.region)
    full_mt = hl.read_matrix_table(args.mt)

    monomorphic_checkpoint = args.output_name + "_non_monomorphic_variants.mt"
    if (not utils.check_exists(monomorphic_checkpoint)) or args.force:
        non_mono_mt = vq.remove_monomorphic(full_mt)
        non_mono_mt = non_mono_mt.checkpoint(monomorphic_checkpoint, overwrite=True)

    else:
        logging.info(f"Detected file with monomorphic variants filtered out: {monomorphic_checkpoint}. Loading this file.")
        non_mono_mt = hl.read_matrix_table(monomorphic_checkpoint)

    start_count = non_mono_mt.count_rows()
    logging.info(f"Number of remaining variants after removing monomorphic variants: {start_count}")

    #########################################################################
    ## Filter failing samples and entries to calculate case/control counts ##
    #########################################################################
    # So not to count 'bad' genotypes/samples in calculating case/control counts for variants.
    filter_checkpoint = args.output_name + "_bad_entries_variants_filtered.mt"

    if (not utils.check_exists(filter_checkpoint)) or args.force:
        mt_filtered = sq.filter_failing(
            non_mono_mt, filter_checkpoint, prefix='final', entries=True, variants=False, samples=True,
            unfilter_entries=False,
            pheno_qc=True, min_dp=10, min_gq=20, max_het_ref_reads=0.8,
            min_het_ref_reads=0.2, min_hom_ref_ref_reads=0.9,
            max_hom_alt_ref_reads=0.1, force=args.force
        )

        mt_filtered = mt_filtered.checkpoint(filter_checkpoint, overwrite=True)
    else:
        logging.info(f"Detected file with bad variants and samples filtered out: {filter_checkpoint}. Loading this file.")
        mt_filtered = hl.read_matrix_table(filter_checkpoint)

    ##################################################
    # Annotate matrix table with case/control counts #
    ##################################################
    # Annotate mt with het/homvar case/control counts
    case_annot_checkpoint = args.output_name + "_carriers_annotated.mt"

    if (not utils.check_exists(case_annot_checkpoint)) or args.force:
        # Do counts with dataset with bad variants and genotypes filtered out
        mt_filt_annot, annotations_to_transfer = cv.count_case_control_carriers(
            mt_filtered, case_annot_checkpoint, pheno_col=args.pheno_col, female_col=args.female_col)

        # Transfer the counts to the unfiltered matrix table
        for annotation in annotations_to_transfer:
            mt_case_count = non_mono_mt.annotate_rows(
                **{annotation: mt_filt_annot.rows()[non_mono_mt.row_key][annotation]}
            )

        mt_case_count = mt_case_count.checkpoint(case_annot_checkpoint, overwrite=True)

    else:
        logging.info(f"Detected file with carriers annotated exists: {case_annot_checkpoint}. Loading this file.")
        mt_case_count = hl.read_matrix_table(case_annot_checkpoint)

    control_count = mt_case_count.aggregate_cols(hl.agg.count_where(mt_case_count[args.pheno_col] == False))
    logging.info(f"Number of controls in dataset: {control_count}")
    case_count = mt_case_count.aggregate_cols(hl.agg.count_where(mt_case_count[args.pheno_col] == True))
    logging.info(f"Number of cases in dataset: {case_count}")
    missing = mt_case_count.aggregate_cols(hl.agg.count_where(~hl.is_defined(mt_case_count[args.pheno_col])))
    logging.info(f"Samples missing case/control information: {missing}")

    if missing > 0:
        logging.info(
            f"Warning- samples missing case/control status will be generally ignored in this pipeline.")

    ###############################################################
    ## Annotate variants with gnomad, mpc, CADD, gnomad mismatch ##
    ###############################################################
    var_annot_checkpoint = args.output_name + "_cadd_mpc_gnomad_annot.mt/"

    if (not utils.check_exists(var_annot_checkpoint)) or args.force:
        mt_var_annot = cv.check_variants_annotated(mt_case_count)
    else:
        mt_var_annot = hl.read_matrix_table(var_annot_checkpoint)

    ###########################################################################
    ## Annotate whether variants are sufficiently rare in gnomad or controls ##
    ###########################################################################
    rare_checkpoint = args.output_name +  "_gnomad_control_rarity_annotated_tmp.mt/"

    if (not utils.check_exists(rare_checkpoint)) or args.force:
        mt_pop_annot = cv.annotate_control_rarity(
            mt_var_annot, rare_checkpoint, args.gnomad_population,
            max_allowed_carrier_dominant=args.max_allowed_carrier_dominant,
            max_allowed_homozygotes_recessive=args.max_allowed_homozygotes_recessive,
            gnomad_AF_cutoff_recessive=args.gnomad_AF_cutoff_recessive
        )

        mt_pop_annot = mt_pop_annot.checkpoint(rare_checkpoint, overwrite=True)

    else:
        logging.info(f"Detected file with boolean columns for variant fulfulling population criteria: {rare_checkpoint}. "
                     f"Loading this file.")
        mt_pop_annot = hl.read_matrix_table(rare_checkpoint)

    ############################
    ## Annotate mt with genes ##
    ############################
    genes_annot_checkpoint = args.output_name + "_genes_annotated.mt/"

    if (not utils.check_exists(genes_annot_checkpoint)) or args.force:
        utils.remove_secondary(args.cluster_name)
        mt_genes_annot = cv.annotate_genes(
            mt_pop_annot, genes_annot_checkpoint
        )  # Triggers shuffles?

        mt_genes_annot = mt_genes_annot.checkpoint(genes_annot_checkpoint, overwrite=True)
    else:
        logging.info(f"Detected file with genes annotated exists: {genes_annot_checkpoint}. Loading that.")
        mt_genes_annot = hl.read_matrix_table(genes_annot_checkpoint)

    if args.disease_genes is not None:
        gene_count = mt_genes_annot.aggregate_rows(hl.agg.counter(mt_genes_annot[args.gene_set_name]))
        logging.info(f"Number of variants in given disease gene set: {gene_count}")

    if args.gnomad_gene_metrics is not None:
        pLI_count = mt_genes_annot.aggregate_rows(hl.agg.counter(mt_genes_annot.high_pLI))
        logging.info(f"Number of variants in high pLI genes (>0.9): {pLI_count}")

    ##########################################
    # Run analysis to find putative variants #
    ##########################################
    h.add_preemptibles(args.cluster_name, args.num_preemptibles)
    variants, filt_mt = find_putative_causal_variants(var_mt, args)
    h.remove_preemptibles(args.cluster_name)
    annotate_denovos_genotypes(variants, filt_mt, args)   # Triggers shuffles?

    ###########################
    # Copy logs and shut down #
    ###########################
    logging.info('Pipeline ran successfully! Copying logs and shutting down cluster in 10 minutes.')
    h.copy_logs_output(args.log_dir, log_file=args.log_file, plot_dir=args.log_dir)
