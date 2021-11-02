"""
Script to calculate kinship with King, plot histogram of kinships, and optionally pull duplicates to new file.
"""
import os
import subprocess
import argparse
import matplotlib.pyplot as plt


def recode_bim(filestem):
    """
    Recode bim file, if needed
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


def run_king(plink_data, plink_path, king_path, remake_bed=False):
    """
    Wrapper function to run King relatedness calculations.
    :param args:
    :return:
    """
    ##########################################################
    # Remake bed files if Hail output has mal-formatted them #
    ##########################################################
    if remake_bed:
        if not os.path.exists(plink_data + '_remade.bed'):
            cmd = f"{plink_path} --bfile {plink_data} --make-bed --out {plink_data + '_remade'}"
            subprocess.call(cmd.split())
            print('Remade bed files.')
        else:
            print('Bed files already remade.')

    #####################
    # Calculate kinship #
    #####################
    if not (os.path.exists(f"{plink_data}.kin") | os.path.exists(f"{plink_data}.kin0")):
        print('Calculating kinship coefficients with King.')
        if remake_bed | os.path.exists(plink_data + '_remade.bed'):
            input_name = plink_data + '_remade.bed'
        else:
            input_name = plink_data + '.bed'

        cmd = f"{king_path} -b {input_name} --kinship --prefix {plink_data} --degree 3 --ibs"
        print(f"King command to run: {cmd}")
        subprocess.call(cmd.split())
    else:
        print("Kinships already calculated, detected files.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs King relatedness calculation and finds maximal independent "
                                                 "set of relatives.")
    parser.add_argument("--plink_data", required=True, help="Plink 1.9 fileset (bed/bim/fam) file stem name.")
    parser.add_argument("--recode_bim", action='store_true', help="Recode bim file to make chrX > X?")
    parser.add_argument("--remake_bed", action='store_true', help="Remake bed files, if Hail output is badly formatted?")
    parser.add_argument("--plink_path", required=True, help="Path of plink executable, if necessary to remake bed files.")
    parser.add_argument("--king_path", required=True, help="Path of King executable.")
    parser.add_argument("--pull_duplicates", action='store_true', help="Pull duplicate samples to separate file?")

    args = parser.parse_args()

    #############################
    # Recode bim file if needed #
    #############################
    if args.recode_bim:
        recode_bim(args.plink_data)

    ################################
    # Run king on input plink file #
    ################################
    run_king(args.plink_data, args.plink_path, args.king_path, remake_bed=args.remake_bed)

    ##########################################################
    # Pull kin and kin0 files and plot histogram of kinships #
    ##########################################################
    for kin_file in [args.plink_data + ".kin0", args.plink_data + ".kin"]:
        if os.path.exists(kin_file):
            kinships = []
            with open(kin_file) as in_f:
                first_line = True
                for line in in_f:
                    words = line.strip().split("\t")
                    if first_line:
                        kin_index = words.index('Kinship')
                        first_line = False
                        continue
                    else:
                        kinships.append(words[kin_index])

            in_f.close()

            n, bins, patches = plt.hist(kinships, bins=50)
            plt.title(f'Histogram of {kin_file} kinships')
            plt.ylabel('Frequency')
            plt.xlabel('Kinship coefficient')
            plt.savefig(f'{kin_file.replace(".", "_")}_hist.png')
        else:
            print(f"{kin_file} does not exist.")

    ############################################################
    # Pull duplicates from kin and kin0 files to separate file #
    ############################################################
    if args.pull_duplicates:
        duplicates = open(f"{args.plink_data}_duplicates.txt", "w")
        duplicates.write("\t".join(["FID1", "ID1", "FID2", "ID2", "N_SNP", "HetHet", "IBS0", "Kinship"]))

        for kin_file in [args.plink_data + ".kin0", args.plink_data + ".kin"]:
            if os.path.exists(kin_file):
                with open(kin_file) as in_f:
                    first_line = True
                    for line in in_f:
                        words = line.strip().split("\t")
                        if first_line:
                            kin_index = words.index('Kinship')
                            if "FID" in words:
                                filetype = "kin0"
                            else:
                                filetype = "kin"
                            first_line = False
                            continue
                        else:
                            if float(words[kin_index]) > 0.354:
                                if filetype == "kin":
                                    new_line = [words[0], words[1], words[0], words[2], words[3], words[6], words[7], words[8]]
                                if filetype == "kin0":
                                    new_line = words

                                duplicates.write("\t".join(new_line))
                in_f.close()
            else:
                print(f"{kin_file} does not exist.")

        duplicates.close()








