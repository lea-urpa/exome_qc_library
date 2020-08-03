"""
This script contains functions to look up variants, samples, and genes from matrix tables.
"""

import hail as hl
import variant_annotation as va

# TODO write these functions
# def get_variant_carriers(mt, variant):
# def filter_to_samples(mt, samples):


def calculate_carrier_counts_ids(mt):
    """
    Calculate for a matrix table the number of het and hom alt genotypes per variant, and the carriers for each.

    :param mt: matrix table to annotate
    :param qc_done: Has QC been done on the dataset? Looks for row annotation 
    :return: 
    """
    mt = mt.annotate_rows(hom_alt_gt_count=hl.agg.count_where(mt.GT.is_hom_var()))
    mt = mt.annotate_rows(hom_alt_carriers=hl.agg.filter(mt.GT.is_hom_var(), hl.agg.collect(mt.s)))
    mt = mt.annotate_rows(het_gt_count=hl.agg.count_where(mt.GT.is_het()))
    mt = mt.annotate_rows(het_carriers=hl.agg.filter(mt.GT.is_het(), hl.agg.collect(mt.s)))

    return mt


def get_lof_carriers(mt):
    # Add a row annotation for the number of hom alt and het carriers if the variant is LOF
    mt = mt.annotate_rows(homalt_lof_count=hl.agg.count_where((mt.LOF == True) & (mt.hom_alt_gt_count > 0)))
    mt = mt.annotate_rows(het_lof_count=hl.agg.count_where((mt.LOF == True) & (mt.het_gt_count > 0)))
    # TODO check whether this produces 0 or missing if mt.LOF is false

    # Count if any of these counts is greater than zero
    lof_homalt_count = mt.aggregate_rows(hl.agg.count_where(mt.homalt_lof_count > 0))
    lof_het_count = mt.aggregate_rows(hl.agg.count_where(mt.het_lof_count > 0))

    # If so, filter to variants that are LOF and the counts are greater than zero
    if (lof_homalt_count > 0) or (lof_het_count > 0):
        print("%s homalt and %s het carriers found." % (lof_homalt_count, lof_het_count))
        lofmt = mt.filter_rows((mt.LOF == True) & ((mt.het_lof_count > 0) | (mt.homalt_lof_count > 0)))

        # Then filter to individuals that are hom alt or het for the remaining samples
        lofmt = lofmt.annotate_cols(alt_gt_count=hl.agg.count_where(lofmt.GT.is_hom_var() | lofmt.GT.is_het()))
        lofmt = lofmt.filter_cols(lofmt.alt_gt_count > 0)

        # Then annotate the carriers for each variant
        lofmt = lofmt.annotate_rows(het_carriers=hl.agg.filter(lofmt.GT.is_het(), hl.agg.collect(mt.s)))
        lofmt = lofmt.annotate_rows(hom_alt_carriers=hl.agg.filter(lofmt.GT.is_hom_var(), hl.agg.collect(lofmt.s)))

        # Then export the rows as a table, print and save to file
        carriers_ht = lofmt.rows()
        carriers_ht.show(-1) # TODO test this
    else:
        carriers_ht = None

    return carriers_ht


#def get_missense_carriers(mt):
    # Test for CADD and MPC annotations


    # Add a new row annotation for the number of hom alt and het carriers if the variant is damaging missense

def get_gene_damaging_carriers(mt, gene, report_all_homozygotes=True):

    # Check if gene/LOF/missense variant annotations have been created by running va.annotate_variants
    try:
        mt.gene.describe()
    except Exception as e:
        print('Missing LOF annotation! Annotating now.')
        print(e)
        try:
            mt = va.annotate_variants(mt)
        except Exception as e:
            print('Error! LOF/missense/synonymous annotation failed. Have you run VEP?')
            print(e)
            return

    # Filter MT to gene of interest
    genemt = mt.filter_rows(mt.gene.contains(gene))

    if genemt.count_rows() == 0:
        print('Gene %s not in dataset. Did you spell it right?' % gene)
        return

    # Calculate carrier counts for individuals in this gene
    genemt = calculate_carrier_counts_ids(genemt)

    # Report damaging variants + carriers
    print('Checking for LOF variants in gene %s' % gene)
    lof_homalt_count, lof_het_count, carriers_ht = get_lof_carriers(genemt)

    print('Saving this table to het_carriers.txt')
    carriers_ht.export("%s_het_carriers.txt" % gene)

    print('Checking for damaging missense variants in gene %s' % gene)

