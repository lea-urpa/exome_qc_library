"""
Script to calculate kinship with King, and then find the maximum set of unrelated individuals.

Author: Lea Urpa, August 2020
"""
import os
import subprocess
import argparse
import networkx as nx


def recode_bim(filestem):
    """
    Takes a Plink filestem and recodes 'chr1' style chromosome IDs to '1' style chromsome IDs.

    :param filestem:
    :return:
    """
    print(f'Recoding bim file: {filestem}.bim')
    old_bim = filestem + "_old_chrom_names.bim"
    new_bim = filestem + ".bim"
    os.rename(new_bim, old_bim)

    new_bim = open(new_bim, "w")
    with open(old_bim) as old:
        for line in old:
            words = line.strip().split('\t')
            new_chrom = words[1].replace("chr", "")
            words[0] = new_chrom

            new_bim.write("\t".join(words) + "\n")


def run_king(args):
    """
    Wrapper function to run King relatedness calculations.
    :param args:
    :return:
    """
    ##########################################################
    # Remake bed files if Hail output has mal-formatted them #
    ##########################################################
    if args.remake_bed:
        if not os.path.exists(args.plink_data + '_remade.bed'):
            cmd = f"{args.plink_path} --bfile {args.plink_data} --make-bed --out {args.plink_data + '_remade'}"
            subprocess.call(cmd.split())
            print('Remade bed files.')
        else:
            print('Bed files already remade.')

    #####################
    # Calculate kinship #
    #####################
    print('Calculating kinship coefficients with King.')
    if args.remake_bed | os.path.exists(args.plink_data + '_remade.bed'):
        input_name = args.plink_data + '_remade.bed'
    else:
        input_name = args.plink_data + '.bed'

    cmd = f"{args.king_path} -b {input_name} --kinship --prefix {args.plink_data} --degree 3"
    print(f"King command to run: {cmd}")
    subprocess.call(cmd.split())


def create_edgelist(filestem):
    """
    Takes king output and creates list of related couples
    :param filestem: filestem of .kin/.kin0 output
    :return:
    """
    print("Formatting king output as edgelist.")
    edgelist = open(filestem + ".edgelist", "w")

    #####################################################
    # Get pairs of related individuals from .kin output #
    #####################################################
    if os.path.isfile(filestem + ".kin"):
        with open(filestem+ ".kin") as kin_in:
            for line in kin_in:
                words = line.strip().split('\t')
                if words[0] != "FID":
                    id1 = words[2]
                    id2 = words[3]

                    edgelist.write("\t".join([id1, id2]) + "\n")
    kin_in.close()

    ######################################################
    # Get pairs of related individuals from .kin0 output #
    ######################################################
    if os.path.isfile(filestem + ".kin0"):
        with open(filestem + ".kin0") as kin0_in:
            for line in kin0_in:
                words = line.strip().split("\t")
                if words[0] != "FID1":
                    id1 = words[2]
                    id2 = words[4]

                    edgelist.write("\t".join([id1, id2]) + "\n")

    kin0_in.close()
    edgelist.close()


def create_case_list(filestem):
    """
    Given a plink filestem, extracts individuals that have affectation status == 2 (case) and collects the sample IDs
    to a list.
    :param filestem: plink filestem name
    :return: returns list of cases
    """
    print(f"Opening {filestem}.fam to get cases to preferentially keep in finding maximally unrelated individuals.")

    cases = []
    with open(filestem + ".fam") as in_fam:
        for line in in_fam:
            words = line.strip().split("\t")
            if words[5] == "2":
                cases.append(words[1])

    in_fam.close()
    print(f"Number of cases: {len(cases)}")

    return cases


def nx_algorithm(g, cases):
    """
    nx native based method for filtering cases. In each subgraph it maximizes the independent cases (read: disease
    affected) first and then proceeds with the rest of the subgraph (by including these cases as required nodes in the
    subsequent maximally independent graph)

    :param g: nx graph
    :param cases: list of cases
    :return: list of related nodes that are discarded
    """
    # Instantiate list of unrelated notes
    unrelated_nodes = []

    ##########################################
    # Loop through subgraphs in larger graph #
    ##########################################
    for subgraph in nx.connected_component_subgraphs(g, copy=True):
        subcases = filter_graph_cases(subgraph, cases)   # list of local cases
        if len(subcases) > 0:
            # graph induced by local cases
            case_graph = subgraph.subgraph(cases)
            # get maximal set of cases in case graph
            unrelated_subcases = nx.maximal_independent_set(case_graph)
            # remove cases that are not in the set of independent cases
            related_subcases = list(set(cases) - set(unrelated_subcases))
            subgraph.remove_nodes_from(related_subcases)
            # get maximal set of nodes in subgraph, giving unrelated cases to keep
            unrelated_nodes += nx.maximal_independent_set(subgraph, unrelated_subcases)
        else:
            unrelated_nodes += nx.maximal_independent_set(subgraph)

    ####################
    # result summaries #
    ####################
    related_nodes = list(set(g.nodes()) - set(unrelated_nodes))
    unrelated_cases = filter_node_cases(unrelated_nodes, cases)

    #####################################################################################
    # sanity_check: makes sure that the unrelated cases/nodes are in fact not connected #
    #####################################################################################
    sanity_check(g, unrelated_nodes)
    sanity_check(g, unrelated_cases)

    return related_nodes, len(unrelated_cases), len(unrelated_nodes)


def filter_graph_cases(graph, cases):
    return [elem for elem in graph.nodes() if elem in cases]


def filter_node_cases(nodes, cases):
    return [elem for elem in nodes if elem in cases]


def sanity_check(g, nodes):
    """
    Given a list of nodes it makes sure that the algorithms are working properly.
    That is, that the subgraph induced by the remaining nodes does not contain edges.
    :param g: network x graph
    :param nodes: notes of a networkx graph
    :return:
    """
    assert g.subgraph(nodes).number_of_edges() == 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs King relatedness calculation and finds maximal independent "
                                                 "set of relatives.")
    parser.add_argument("--plink_data", required=True, help="Plink 1.9 fileset (bed/bim/fam) file stem name.")
    parser.add_argument("--skip_king", action='store_true', help="Skip calculating kinship coefficients.")
    parser.add_argument("--king_path", help="Path of plink executable.")
    parser.add_argument("--plink_path", help="Path of plink executable, if remaking bed needed.")
    parser.add_argument("--remake_bed", action='store_true', help="Remake bed files?")
    parser.add_argument("--recode_bim", action='store_true', help="Recode bim file chromosomes?")

    arguments = parser.parse_args()

    ################################
    # Recode bim file if necessary #
    ################################
    if arguments.recode_bim:
        recode_bim(arguments.plink_data)

    ################################
    # Run king on input plink file #
    ################################
    if not arguments.skip_king:
        if arguments.king_path is None:
            print("Error! --king_path must be given if not skipping king calculations.")
        run_king(arguments)
    else:
        if not (os.path.isfile(arguments.plink_data + ".kin") | (os.path.isfile(arguments.plink_data + ".kin0"))):
            print(f"Error! File {arguments.plink_data + '.kin'} or {arguments.plink_data + '.kin0'} must exist if "
                  f"you are skipping kinship coefficient calculation.")
            exit()

    #########################################################
    # Format kin output to edgelist, get list of cases only #
    #########################################################
    create_edgelist(arguments.plink_data)
    case_list = create_case_list(arguments.plink_data)

    ##################################################
    # Input as graph and get maximal independent set #
    ##################################################
    print(f'Input file for related cases: {arguments.plink_data + ".edgelist"}')
    related_ind_g = nx.read_edgelist(arguments.plink_data + ".edgelist")

    print('Calculating maximal independent set.')
    related_to_remove, num_ind_cases, num_ind_nodes = nx_algorithm(related_ind_g, case_list)
    print('# of unrelated cases given king input: ' + str(num_ind_cases))

    print('Writing file with related individuals to remove from dataset:')
    with open('related_to_remove.txt', 'w') as f:
        for item in related_to_remove:
            f.write("%s\n" % item)
