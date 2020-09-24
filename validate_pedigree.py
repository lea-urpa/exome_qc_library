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
from collections import OrderedDict
import utils


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


def parse_kin_errors(kin_file, cohorts, args):
    """
    Read in .kin file and report errors in readable way
    :param kin_file: .kin output from King
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
        delim = utils.check_delimiter(kin_file)
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


def parse_kin0_errors(kin0_file, cohorts, args):
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

        delim = utils.check_delimiter(kin0_file)
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


def count_families_kin(kin_file, cohorts, kinship_cutoff):
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

        delim = utils.check_delimiter(kin_file)
        with open(kin_file) as kin_in:
            for line in kin_in:
                if not ((cohort == "all") | (cohort in line)):
                    continue

                fid, id1, id2, n_snp, z0, phi, hethet, ibs0, kinship, error = line.strip().split(delim)
                if (fid != "FID") and (fid != "0") and (error == "0"):
                    edgelist.write("\t".join([id1, id2]) + "\n")

        kin_in.close()

        ###################################################################
        # Read in the edgelist and count the number of independent graphs #
        ###################################################################
        g = nx.read_edgelist(edgelist_name)

        num_families = nx.number_connected_components(g)

        families_count = {}
        for c in nx.connected_components(g):
            node_count = len(c)
            if node_count not in families_count.keys():
                families_count[node_count] = 1
            else:
                families_count[node_count] += 1
            # If the node count is too high, write names to file
            if node_count > 50:
                print("Problematic family detected- finding individuals connected to too many others and writing "
                      "names to file.")
                problem_graph = g.subgraph(c)
                most_connected = max(dict(problem_graph.degre()).items(), key=lambda x: x[1])

                filename = f"family_{node_count}_members"
                while os.path.isfile(filename):
                    filename = filename + "-1"
                filename = filename + ".txt"
                out_f = open(filename, "w")
                out_f.write("ID\tnum_connections\n")

                connections = dict(problem_graph.degree())
                connections_ordered = OrderedDict(sorted(connections.items(), key=lambda kv: kv[1], reverse=True))

                for id, num_connections in connections_ordered.items():
                    out_f.write(f"{id}\t{num_connections}\n")

        families_count = OrderedDict(sorted(families_count.items()))
        for family_size, count in families_count.items():
            print(f"{count} families with {family_size} individuals.")


def count_trios(kin_file, fam_file,args):
    pass
    ############################################################################
    # Read in kin + fam file and count the number of complete, validated trios #
    ############################################################################

    ####################################################
    # Read in fam file and collect dictionary of trios #
    ####################################################



# Create edgelist from kin file without errors and kin0 file and report # of families, and size of families (subtracting

# individuals who are related to too many other individuasl > likely mixups)

