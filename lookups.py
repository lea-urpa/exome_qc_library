"""
Contains functions to look up variants, samples, and genes from matrix tables.

Author: Lea Urpa, August 2020
"""


def calculate_carrier_counts_ids(mt):
    """
    Calculate for a matrix table the number of het and hom alt genotypes per variant, and the carriers for each.

    :param mt: matrix table to annotate
    :return: 
    """
    mt = mt.annotate_rows(hom_ref_gt_count=hl.agg.count_where(mt.GT.is_hom_ref() & hl.is_defined(mt.GT)))
    mt = mt.annotate_rows(het_gt_count=hl.agg.count_where(mt.GT.is_het() & hl.is_defined(mt.GT)))
    mt = mt.annotate_rows(het_carriers=hl.agg.filter(mt.GT.is_het() & hl.is_defined(mt.GT), hl.agg.collect(mt.s)))
    mt = mt.annotate_rows(hom_alt_gt_count=hl.agg.count_where(mt.GT.is_hom_var() & hl.is_defined(mt.GT)))
    mt = mt.annotate_rows(hom_alt_carriers=hl.agg.filter(mt.GT.is_hom_var() & hl.is_defined(mt.GT),
                                                         hl.agg.collect(mt.s)))

    return mt


def get_lof_carriers(mt, args):

    ###############################################################
    # Count the number of LOF homalt or het genotypes in the gene #
    ###############################################################
    lof_homalt_ct = mt.aggregate_rows(
        hl.agg.count_where(hl.is_defined(mt.LOF) & (mt.LOF == True) & (mt.hom_alt_gt_count > 0)))
    lof_het_ct = mt.aggregate_rows(
        hl.agg.count_where(hl.is_defined(mt.LOF) & (mt.LOF == True) & (mt.het_gt_count > 0)))

    ################################################################
    # If any, filter mt to those variants and export to Hail table #
    ################################################################
    if (lof_homalt_ct > 0) or (lof_het_ct > 0):
        print(f"{lof_homalt_ct} homalt and {lof_het_ct} het carriers found.")

        lofmt = mt.filter_rows((mt.LOF == True) & ((mt.het_gt_count > 0) | (mt.hom_alt_gt_count > 0)), keep=True)

        # Then export the rows as a table, print and save to file
        carriers_ht = lofmt.rows()
        carriers_ht.export(f"{args.output_stem}_{gene}_LOF_carriers.txt")

    else:
        print(f"No LOF variant carriers found in gene {gene}")


def get_variant_carriers(mt, args):

    #################################################
    # Find carriers for variant if it is in dataset #
    #################################################
    rowct = mt.count_rows()
    if rowct == 0:
        print(f"No variants in list in dataset.")
        exit()

    ############################
    # Calculate carrier counts #
    ############################
    mt = calculate_carrier_counts_ids(mt)
    if args.pheno_col is not None:
        mt = count_case_control_carriers(mt, args)

    #####################################
    # Pull variants to rows, checkpoint #
    #####################################
    rows = mt.rows()
    rows = rows.checkpoint(args.output_stem + "_rows_tmp.mt", overwrite=True)
    rows.write(f"{args.output_stem}_carriers.ht/", overwrite=True)

    cols = mt.cols()

    rows = rows.annotate(popmax_population=rows.gnomad_popmax[rows.gnomad_popmax_index_dict['gnomad']].pop)
    rows = rows.annotate(gnomad_fin_freq=rows.gnomad_freq[rows.gnomad_freq_index_dict['gnomad_fin']],
                         gnomad_controls_fin_freq=rows.gnomad_freq[rows.gnomad_freq_index_dict['controls_fin']],
                         gnomad_nonneuro_fin_freq=rows.gnomad_freq[rows.gnomad_freq_index_dict['non_neuro_fin']],
                         gnomad_popmax_gnomad=rows.gnomad_popmax[rows.gnomad_popmax_index_dict['gnomad']],
                         gnomad_popmax_controls=rows.gnomad_popmax[rows.gnomad_popmax_index_dict['controls']],
                         gnomad_popmax_nonneuro=rows.gnomad_popmax[rows.gnomad_popmax_index_dict['non_neuro']])

    rows = rows.drop('gnomad_freq', 'gnomad_popmax', 'vep')

    #############################################################################
    # Drop hom alt carriers and explode by het carriers, flatten, write to file #
    #############################################################################
    het_carriers = rows.drop("hom_alt_carriers")
    het_carriers = het_carriers.explode(het_carriers.het_carriers)
    het_carriers = het_carriers.flatten()
    het_carriers.export(f"{args.output_stem}_het_carriers.txt")

    #############################################################################
    # Drop het carriers and explode by hom alt carriers, flatten, write to file #
    #############################################################################
    homalt_carriers = rows.drop("het_carriers")
    homalt_carriers = homalt_carriers.explode(homalt_carriers.hom_alt_carriers)
    homalt_carriers = homalt_carriers.flatten()
    homalt_carriers.export(f"{args.output_stem}_homalt_carriers.txt")

    ########################################################
    # Export columns for the individuals in the hail table #
    ########################################################
    cols = cols.flatten()
    cols.export(f"{args.output_stem}_carrier_info.txt")


if __name__ == "__main__":
    ######################################
    # Initialize Hail and import scripts #
    ######################################
    import os
    import argparse
    import hail as hl
    hl.init()

    parser = argparse.ArgumentParser(description="Lookups for damaging variants in genes from exome sequencing data.")
    parser.add_argument("-mt", help="Matrix table to search for damaging variants in genes and/or specific variants.")
    parser.add_argument("--pheno_col", type=str, help="Col annotation giving true/false for case phenotype")
    parser.add_argument("--female_col", type=str, help="Col annotation giving true/false for female status.")
    parser.add_argument("--output_name", required=True, type=str, help="Output name stem for results.")
    parser.add_argument("--output_dir", required=True, type=str, help="Output directory for results.")
    parser.add_argument("--genes", type=str,
                        help="comma-separated genes or single gene in which to search for damaging variants.")
    parser.add_argument("--variants", type=str, help="comma-separated variants or single variant to report carriers.")
    parser.add_argument("--variant_list", type=str,
                        help="file containing list of variants, one per line, to report carriers.")
    parser.add_argument("--variant_sep", type=str, default=":", help="Separator between chrom/pos/ref/alt")
    parser.add_argument("--scripts_dir", required=True, help="Directory for exome qc library scripts")
    parser.add_argument("--reference_genome", default="GRCh38", choices=["GRCh37",  "GRCh38"])
    parser.add_argument("--cadd_ht", help="Hail table with CADD variant information")
    parser.add_argument("--mpc_ht", help="Hail table with MPC varian information")
    parser.add_argument("--gnomad_ht", help="Hail table with gnomad variant information")
    parser.add_argument("--gnomad_mismatch_ht", help="Hail table with variant gnomad mismatch information")

    args = parser.parse_args()

    if (args.genes is None) and (args.variants is None) and (args.variant_list is None):
        print("Error! One of --genes, --variants, or --variant_list must be given!")
        exit()

    if (args.pheno_col is not None) and (args.female_col is None):
        print("Error! if giving --pheno_col, --female_col must also be given")

    args.output_stem = os.path.join(args.output_dir, args.output_name)

    ##########################
    # Import python scripts  #
    ##########################
    scripts = ["variant_annotation.py", "find_putative_causal_variants.py"]
    for script in scripts:
        hl.spark_context().addPyFile(os.path.join(args.scripts_dir, script))

    import variant_annotation as va
    from find_putative_causal_variants import count_case_control_carriers

    ########################
    # Load in matrix table #
    ########################
    fullmt = hl.read_matrix_table(args.mt)

    ####################################
    # Check if variant annotation done #
    ####################################
    try:
        fullmt.gene.describe()
    except Exception as e:
        print('Missing LOF annotation! Annotating now.')
        print(e)
        try:
            fullmt = va.annotate_variants(fullmt)
        except Exception as e:
            print('Error! LOF/missense/synonymous annotation failed. Have you run VEP?')
            print(e)
            exit()

    if args.genes is not None:
        #################################
        # Filter MT to gene of interest #
        #################################
        gene_list = args.genes.strip().split(",")

        for gene in gene_list:
            genemt = fullmt.filter_rows(fullmt.gene.contains(gene))
            genemt = genemt.checkpoint(args.output_stem + "_genemt_tmp.mt", overwrite=True)

            if genemt.count_rows() == 0:
                print(f'Gene {gene} not in dataset. Did you spell it right?')
                continue

            #########################################################
            # Calculate carrier counts for individuals in this gene #
            #########################################################
            genemt = calculate_carrier_counts_ids(genemt)

            #######################################
            # Report damaging variants + carriers #
            #######################################
            print(f'Checking for LOF variants in gene {gene}')
            get_lof_carriers(genemt, args)

    if (args.variants is not None) or (args.variant_list is not None):
        ##############################################################
        # Import list of variants, from command line or list or both #
        ##############################################################
        variants = []
        if args.variant_list is not None:
            variant_table = hl.import_table(args.variant_list, no_header=True)
            variants = variant_table.f0.take(variant_table.count())
            variants = [x.split(args.variant_sep) for x in variants]

        if args.variants is not None:
            args_vars = args.variants.strip().split(",")
            args_vars = [x.split(args.variant_sep) for x in args_vars]
            variants.extend(args_vars)
        print(f"Number of variants to find carriers for: {len(variants)}")

        ##################################################################################
        # Get array of loci in Hail format, do an intermediate filter to just those loci #
        ##################################################################################
        if args.reference_genome == "GRCh37":
            hail_loci = hl.array([hl.struct(locus=hl.locus(x[0], int(x[1]), reference_genome="GRCh37"),
                                            alleles=hl.array([x[2], x[3]])) for x in variants])
        elif args.reference_genome == "GRCh38":
            hail_loci = hl.array([hl.struct(locus=hl.locus('chr' + x[0], int(x[1]), reference_genome="GRCh38"),
                                            alleles=hl.array([x[2], x[3]])) for x in variants])
        else:
            print("Incorrect reference genome given! How did we get here?")
            print(args.reference_genome)
            exit()

        variants_mt = fullmt.filter_rows(hail_loci.contains(hl.struct(locus=fullmt.locus, alleles=fullmt.alleles)))

        ######################################################################
        # Annotate intermediate matrix table with annotation files, if given #
        ######################################################################
        if args.cadd_ht is not None:
            print("Anotating variants with CADD info")
            variants_mt = va.annotate_variants_cadd(variants_mt, args.cadd_ht)
        if args.mpc_ht is not None:
            print("Annotating variants with MPC info")
            variants_mt = va.annotate_variants_mpc(variants_mt, args.mpc_ht)
        if args.gnomad_ht is not None:
            print("Annotating variants with gnomad info")
            variants_mt = va.annotate_variants_gnomad(variants_mt, args.gnomad_ht)
        if args.gnomad_mismatch_ht is not None:
            print("Annotating variants with gnomad mismatch info")
            variants_mt = va.annotate_variants_gnomad_mismatch(variants_mt, args.gnomad_mismatch_ht)

        #################################################################
        # Loop through variants and find carriers, export separate file #
        #################################################################
        get_variant_carriers(variants_mt, args)

