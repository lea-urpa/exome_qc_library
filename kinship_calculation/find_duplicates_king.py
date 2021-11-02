"""
Script to calculate kinship with King
"""


import os
import subprocess
import argparse



    :param filestem:
    :return:
    """
    print(f'Recoding bim file: {filestem}.bim')
    old_bim = filestem + "_old_chrom_names.bim"
    new_bim = filestem + ".bim"
    if not os.path.isfile(old_bim):
        os.rename(new_bim, old_bim)

    new_bim = open(new_bim, "w")
    with open(old_bim) as old:
        for line in old:
            words = line.strip().split('\t')
            new_chrom = words[0].replace("chr", "")
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

    cmd = f"{args.king_path} -b {input_name} --kinship --prefix {args.plink_data} --degree 3 --ibs"
    print(f"King command to run: {cmd}")
    subprocess.call(cmd.split())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs King relatedness calculation and finds maximal independent "
                                                 "set of relatives.")
    parser.add_argument("--plink_data", required=True, help="Plink 1.9 fileset (bed/bim/fam) file stem name.")
    parser.add_argument("--kinship_cutoff", default=0.0883,
                        help="Kinship threshold to declare a pair of related individuals, of which one will be removed."
                             " Default 0.0883, approximately the lower cutoff for as second degree relative.")
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