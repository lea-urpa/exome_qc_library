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
import shutil


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
                        if error == "0.5":
                            error_type = "warning"
                        elif error == "1":
                            error_type = "error"
                        else:
                            print(f"Error- unexpected error type: {error}")
                            print(f"line: {line}")
                            exit()

                        line = f"Kinship {error_type} in family {fid} between samples {id1} and {id2}." \
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

        error_report.close()
        in_kin0.close()


def count_families(edgelist_name):
    """
    Given an edgelist (list of pairs of relatives), counts the number of independent graphs (families) and the number
    of nodes in each graph (number of individuals in each family). Also reports unusually large families.
    :param edgelist_name: File name of edgelist ot read
    :return:
    """
    ############################################
    # Read in edgelist, get number of families #
    ############################################
    g = nx.read_edgelist(edgelist_name)

    num_families = nx.number_connected_components(g)
    print(f"Number of families detected from edgelist: {num_families} (file {edgelist_name})")

    # TODO add check for the number of connections each individual has- if it's greater than X, write that in a
    #  report too

    #######################################################
    # Collect dictionary of # of families by size, report #
    #######################################################
    families_count = {}
    for c in nx.connected_components(g):
        node_count = len(c)
        if node_count not in families_count.keys():
            families_count[node_count] = 1
        else:
            families_count[node_count] += 1

        ######################################################
        # If the node count is too high, write names to file #
        ######################################################
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

            out_f.close()

    families_count = OrderedDict(sorted(families_count.items()))
    for family_size, count in families_count.items():
        print(f"{count} families with {family_size} individuals.")


def count_families_kin(kin_file, cohorts):
    """
    Create edgelist from kin file without errors and report # of families, and size of families
    :param kin_file: Kinship output file from King
    :param cohorts: Cohort substrings to subset report to.
    :return:
    """
    if "all" not in cohorts:
        cohorts.append("all")

    for cohort in cohorts:
        if cohort == "all":
            print(f"Looking for # families in all individuals, only in reported, validated families.")
        else:
            print(f"Looking for # of families in cohort {cohort}, only in reported, validated families." )
        ##############################################################################
        # Create an edgelist containing just the expected families, minus the errors #
        ##############################################################################
        filestem = os.path.splitext(kin_file)
        edgelist_name = f"{filestem}_{cohort}.edgelist_kin"
        edgelist = open(edgelist_name, "w")

        delim = utils.check_delimiter(kin_file)
        with open(kin_file) as kin_in:
            for line in kin_in:
                if not ((cohort == "all") or (cohort in line)):
                    continue

                fid, id1, id2, n_snp, z0, phi, hethet, ibs0, kinship, error = line.strip().split(delim)
                if (fid != "FID") and (fid != "0") and (error == "0"):
                    edgelist.write("\t".join([id1, id2]) + "\n")

        kin_in.close()
        edgelist.close()

        ###################################################################
        # Read in the edgelist and count the number of independent graphs #
        ###################################################################
        count_families(edgelist_name)


def count_families_kin0(kin0_file, cohorts):
    # Create edgelist from kin file without errors and kin0 file and report # of families, and size of families
    # (subtracting individuals who are related to too many other individuasl > likely mixups)

    if "all" not in cohorts:
        cohorts.append("all")

    for cohort in cohorts:
        if cohort == "all":
            print(f"Looking for # families in all individuals, including unexpected relationships")
        else:
            print(f"Looking for # of families in cohort {cohort}, including unexpected relationships.")

        #####################################################################
        # Copy old edgelist to new file, append unexpected kinships to that #
        #####################################################################
        filestem = os.path.splitext(kin0_file)
        edgelist_name = f"{filestem}_{cohort}.edgelist.kin"
        new_edgelist_name = edgelist_name + "0"

        shutil.copy(edgelist_name, new_edgelist_name)
        new_edgelist = open(new_edgelist_name, "a")

        delim = utils.check_delimiter(kin0_file)
        with open(kin0_file) as kin0_in:
            for line in kin0_in:
                if not ((cohort == "all") or (cohort in line)):
                    continue

                fid1, id1, fid2, id2, n_snp, hethet, ibs0, kinship = line.strip().split(delim)
                if fid1 != "FID1":
                    new_edgelist.write("\t".join([id1, id2]) + "\n")

        kin0_in.close()
        new_edgelist.close()

        ###############################################################
        # Read in edgelist and count hte number of independent graphs #
        ###############################################################
        count_families(new_edgelist_name)

def count_trios(kin_file, fam_file, cohorts):

    if "all" not in cohorts:
        cohorts.append("all")

    ############################################################################
    # Read in kin + fam file and count the number of complete, validated trios #
    ############################################################################
    for cohort in cohorts:
        if cohort == "all":
            print("Finding trios in whole dataset")
        else:
            print(f"Finding trios in cohort {cohort}")

        ##########################################
        # Read in fam file to get supposed trios #
        ##########################################
        raw_trios = []
        delim1 = utils.check_delimiter(fam_file)
        with open(fam_file) as fam_in:
            for line in fam_in:
                fid, iid, pid, mid, sex, pheno = line.strip().split(delim1)
                if (pid != "0") and (mid != "0"):
                    raw_trios.append({'proband': iid, 'mother': mid, 'father': pid})

        fam_in.close()
        print(f"Number of expected trios from pedigree: {len(raw_trios)} (file {fam_file}")

        #################################################################################
        # Collect sample IDs in kinship file as set, filter raw trios to complete trios #
        #################################################################################
        present_ids = set()
        delim2 = utils.check_delimiter(kin_file)
        with open(kin_file) as kin_in:
            for line in kin_in:
                fid, id1, id2, n_snp, z0, phi, hethet, ibs0, kinship, error = line.strip().split(delim2)

                if fid != "FID":
                    present_ids.add(id1)
                    present_ids.add(id2)

        kin_in.close()

        complete_trios = []
        for trio in raw_trios:
            if (trio['mother'] in present_ids) and (trio['father'] in present_ids) and (trio['proband'] in present_ids):
                trio['m_valid'] = False
                trio['p_valid'] = False
                complete_trios.append(trio)

        print(f"Number of complete trios from pedigree: {len(raw_trios)} (sample in kinship data)")

        ###################################################################################
        # Loop through kinship file and see if matching PO in trios exists, annotate trio #
        ###################################################################################
        # This might take a long time, but it should work
        with open(kin_file) as kin_in:
            for line in kin_in:
                fid, id1, id2, n_snp, z0, phi, hethet, ibs0, kinship, error = line.strip().split(delim2)

                for trio in complete_trios:
                    if ((id1 == trio['proband']) and (id2 == trio['mother'])) or \
                       ((id2 == trio['proband']) and (id1 == trio['mother'])):
                        relationship = infer_relationship(ibs0, kinship)
                        if relationship == "parent-offspring":
                            trio['m_valid'] == True

                    if ((id1 == trio['proband']) and (id2 == trio['father'])) or \
                       ((id2 == trio['proband']) and (id1 == trio['father'])):
                        relationship = infer_relationship(ibs0, kinship)
                        if relationship == "parent-offspring":
                            trio['p_valid'] == True

        valid_trios = []
        for trio in complete_trios:
            if (trio['m_valid'] == True) and (trio['p_valid'] == True):
                valid_trios.append(trio)

        print(f"Number of valid trios: {len(valid_trios)}")


# Arguments and file main goes here
