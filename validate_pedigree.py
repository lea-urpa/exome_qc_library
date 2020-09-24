"""
Scripts to take King genetic relatedness information + the pedigree/fam file that went into King, and creates
a report of the errors and the number of families/sizes of families in the dataset.
"""
# Needed args
# cohort_subset_string
# output_stem (containing output name)
# cohorts = args.cohort_subset_string.strip().split(",")

import time
import os
import networkx as nx
from calculate_kinship_king import create_edgelist


def infer_relationship(ibs0, kinship):
    """
    Given ibs0 and kinship values from King output, infers relationship up to 3rd degree.
    :param ibs0:
    :param kinship:
    :return: returns string giving inferred relationship
    """
    if kinship > 0.354:
        actual = "MZ twin or duplicate"
    elif (kinship <= 0.354) and (kinship > 0.177) and (ibs0 < 0.1):
        actual = "parent-offspring"
    elif (kinship <= 0.354) and (kinship > 0.177) and (ibs0 >= 0.1):
        actual = "full sibling"
    elif (kinship <= 0.177) and (kinship > 0.0884):
        actual = "half sibling, avuncular, or grandparent-grandchild"
    elif (kinship <= 0.0884) and (kinship > 0.0442):
        actual = "first cousin or half avuncular"
    else:
        actual = "unrelated"

    return actual


def parse_kin_errors(kin_file, delim, cohorts, args):
    """
    Read in .kin file and report errors in readable way
    :param kin_file: .kin output from King
    :param delim: delimiter for King output
    :param cohorts: array of strings with cohort subset strings to subset report to.
    :param args: output stem to write to
    :return:
    """
    if "all" not in cohorts:
        cohorts.append("all")

    ########################
    # Loop through cohorts #
    ########################
    for cohort in cohorts:
        error_report = open(f"{args.output_stem}_{cohort}_kin_errors.txt", "w")
        error_report.write(f"#Input kinship file: {kin_file}\n")
        error_report.write(f"#Date and time: {time.asctime()}\n")
        kinship_errors = 0
        with open(kin_file) as in_kin:
            for line in in_kin:
                # Create separate report for all errors in kinship file, or just IDs containing cohort of interest
                if not ((cohort == "all") | (cohort in line)):
                    continue

                fid, id1, id2, n_snp, z0, phi, hethet, ibs0, kinship, error = line.strip().split(delim)
                if fid != "FID":
                    z0 = float(z0)
                    phi = float(phi)
                    ibs0 = float(ibs0)
                    kinship = float(kinship)
                    if error != "0":  # marked as error in kinship
                        kinship_errors += 1

                        ################################################
                        # Infer expected relationships from z0 and phi #
                        ################################################
                        if (z0 == 0) and (phi == 0.5):
                            expected = "MZ twin or duplicate"
                        elif (z0 == 0) and (phi == 0.25):
                            expected = "parent-offspring"
                        elif (z0 == 0.25) and (phi == 0.25):
                            expected = "full sibling"
                        elif (z0 == 0.5) and (phi == 0.125):
                            expected = "half sibling, avuncular, or grandparent-grandchild"
                        elif (z0 == 0.75) and (phi == 0.0625):
                            expected = "first cousin or half-avuncular"
                        elif (z0 == 1) and (phi == 0):
                            expected = "unrelated"
                        else:
                            print(f"Error! Unexpected expected kinship/ibs0 values. z0: {z0} phi: {phi}")
                            exit()

                        ####################################################
                        # Infer actual relationships from ibs0 and kinship #
                        ####################################################
                        actual = infer_relationship(ibs0, kinship)

                        line = f"Kinship error in family {fid} between samples {id1} and {id2}." \
                            f"Expected relationship: {expected} Kinship supports: {actual}. "\
                            f"IBS0: {ibs0}. Kinship: {kinship} \n"
                        error_report.write(line)

        error_report.close()
        in_kin.close()


def parse_kin0_errors(kin0_file, delim, cohorts, args):
    """
    Read in .kin0 file and report errors in readable way
    :param kin0_file: .kin0 file output from King
    :param delim: delimiter for .kin0 file
    :param cohorts: cohorts to subset reports to
    :param args: arguments
    :return:
    """
    if "all" not in cohorts:
        cohorts.append("all")

    for cohort in cohorts:
        error_report = open(f"{args.output_stem}_{cohort}_kin0_errors.txt", "w")
        error_report.write(f"#Input kinship file: {kin0_file}\n")
        error_report.write(f"#Date and time: {time.asctime()}\n")

        with open(kin0_file) as in_kin0:
            for line in in_kin0:
                fid1, id1, fid2, id2, n_snp, hethet, ibs0, kinship = line.strip().split(delim)

                ####################################################
                # Infer actual relationships from ibs0 and kinship #
                ####################################################
                actual = infer_relationship(ibs0, kinship)

                line = f"Unexpected kinship between {id1} ({fid1}) and {id2} ({fid2}). Expected to be unrelated; kinship " \
                    f"supports {actual}. IBS0: {ibs0}. Kinship: {kinship}"

                error_report.write(line)


def count_families_kin(kin_file, delim, cohorts, args):
    # Create edgelist from kin file without errors and report # of families, and size of families
    if "all" not in cohorts:
        cohorts.append("all")

    for cohort in cohorts:
        ##############################################################################
        # Create an edgelist containing just the expected families, minus the errors #
        ##############################################################################
        filestem = os.path.splitext(kin_file)
        edgelist_name = f"{filestem}_{cohort}.edgelist_kin"
        edgelist = open(edgelist_name, "w")

        with open(kin_file) as kin_in:
            for line in kin_in:
                if not ((cohort == "all") | (cohort in line)):
                    continue

                fid, id1, id2, n_snp, z0, phi, hethet, ibs0, kinship, error = line.strip().split(delim)
                if fid != "FID":
                    if error == "0":
                        edgelist.write("\t".join([id1, id2]) + "\n")

        ###################################################################
        # Read in the edgelist and count the number of independent graphs #
        ###################################################################
        g = nx.read_edgelist(edgelist_name)





# Create edgelist from kin file without errors and kin0 file and report # of families, and size of families (subtracting

# individuals who are related to too many other individuasl > likely mixups)

