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
        carriers_ht.export(f"{args.output_stem}_{gene}LOF_carriers.txt")

    else:
        print(f"No LOF variant carriers found in gene {gene}")


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
    parser.add_argument("--output_name", required=True, type=str, help="Output name stem for results.")
    parser.add_argument("--output_dir", required=True, type=str, help="Output directory for results.")
    parser.add_argument("--genes", required=True, type=str,
                        help="comma-separated genes or single gene in which to search for damaging variants.")
    parser.add_argument("--scripts_dir", required=True, help="Directory for exome qc library scripts")
    parser.add_argument("--reference_genome", default="GRCh38", choices=["GRCh37",  "GRCh38"])

    args = parser.parse_args()

    args.output_stem = os.path.join(args.output_dir, args.output_name)

    ##########################
    # Import python scripts  #
    ##########################
    scripts = ["variant_annotation.py"]
    for script in scripts:
        hl.spark_context().addPyFile(os.path.join(args.scripts_dir, script))

    import variant_annotation as va

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

    #################################
    # Filter MT to gene of interest #
    #################################
    gene_list = args.genes.strip().split(",")

    for gene in gene_list:
        genemt = fullmt.filter_rows(fullmt.gene.contains(gene))

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

