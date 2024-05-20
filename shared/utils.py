"""
utility functions
"""
import os
from functools import partial
from tempfile import NamedTemporaryFile
import gzip
import csv
import subprocess
import shlex
import time
import logging
import hail as hl


def cluster_region_apriori(cluster_name):
    """
    Given the name of a cluster, pulls the cluster region from description output. Takes a few seconds for each check.
    :param cluster_name: name of dataproc cluster
    :return: cluster region name
    """
    # Get cluster region
    cluster_output = None
    for region in ["europe-west1", "europe-north1", "europe-west2", "europe-west3", "europe-west6", "europe-central2",
                   "us-central1", "us-east1", "us-east4", "us-west1", "us-west2", "us-west3", "us-west4"]:

        cluster_describe = f"gcloud dataproc clusters describe {cluster_name} --region {region}"

        try:
            cluster_output = subprocess.check_output(shlex.split(cluster_describe))
        except:
            pass

    if cluster_output == None:
        print("Error! Cluster not found. Currently only europe and us regions supported for this function.")

    cluster_clean = str(cluster_output).replace("\\t", "").split("\\n")
    cluster_region = [x.strip() for x in cluster_clean if "goog-dataproc-location" in x][0].split(":")[1].strip()

    return(cluster_region)


def check_regions(cluster_region, file_url):
    """
    Checks that cluster region and region of the file inputs or outputs are the continent, or stops the pipeline.
    :param cluster_region: string input, cluster region
    :param file_url: string input, file to read in or out
    :return:
    """
    # Get bucket region
    if not file_url.startswith("gs://"):
        print("Error! File url input does not start with gs://. Is this a google bucket url?")
        exit()
    bucket = "gs://" + file_url.split("/")[2]
    bucket_describe = f"gsutil ls -L -b {bucket}"

    try:
        bucket_output = subprocess.check_output(shlex.split(bucket_describe))
    except Exception as e:
        print(f"Error! File {file_url} does not appear to exist or is not in Google Cloud Storage.")
        print(e)

    bucket_clean = str(bucket_output).replace("\\t", "").split("\\n")
    bucket_region = [x for x in bucket_clean if "location constraint" in x.lower()][0].split(":")[1].split("-")[0].lower()

    cluster_continent = cluster_region.split("-")[0].lower()

    if cluster_continent != bucket_region:
        print("Error! Bucket region is in a different continent than dataproc cluster. Stoppping pipeline!")
        exit()


def load_vcfs(vcf_files, data_dir, out_dir, combined_mt_fn, force=False, test=False, chr_prefix=False,
              reference_genome="GRCh38", force_bgz=False, force_load=False, call_fields="PGT", save_row_annots=False,
              chrom_split=False, sample_split=False):
    """
    Given a list of VCF files, imports and combines them into a single matrix table.
    :param vcf_files: list of strings (vcf files names)
    :param data_dir: string, directory (in the cloud or local) in which the VCF files are located.
    :param out_dir: string, location to save imported data
    :param combined_mt_fn: combined matrix table file name
    :param force: Force re-import even if files exist?
    :param test: Filter to chr22 for each file when importing?
    :param chr_prefix: Recodes chromosome prefix (e.g. Chr11) while loading data if mismatching with reference genome.
    :param reference_genome: either GRCh38 or GRCh37, exactly.
    :param force_bgz: Force block gzip load
    :param force_load: Forces loading of standard (non-block) gzipped file. Takes longer!
    :param call_fields: Call field for VCF file
    :param save_row_annots: Save variant (row) annotations from multiple input VCFs?
    :param chrom_split: Is the dataset split by chromosome? Must have exactly same samples in each file.
    :param sample_split: Are all your VCFs single-sample VCFS? If so, takes wildcard to get list of matching VCF files.
    """
    logging.debug(
        f"Load vcfs parameters:\n"
        f"vcf_files: {vcf_files}\n"
        f"data_dir: {data_dir}\n"
        f"out_dir: {out_dir}\n"
        f"force reloading vcfs?: {force}\n"
        f"test: {test}\n"
        f"Fix chromosome prefixes? {chr_prefix}\n"
        f"reference genome: {reference_genome}\n"
        f"force block gzip loading? {force_bgz}\n"
        f"force load non-block-gzipped (standard gzipped) file?: {force_load}\n"
        f"call_fields: {call_fields}\n"
        f"save row annotations? {save_row_annots}\n"
        f"\n"
    )

    # Get combined mt output name
    if test:
        test_str = "_test"
    else:
        test_str = ""

    logging.info('Importing genotype vcfs.')
    matrix_tables = []

    ################################################################
    # Deal with mismatching chromosome codes with reference genome #
    ################################################################
    if (chr_prefix is True) and (reference_genome == "GRCh37"):
        recode = {f"chr{i}": f"{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}

    elif (chr_prefix is False) and (reference_genome == "GRCh38"):
        recode = {f"{i}": f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
    else:
        recode = None

    ############################################################
    # If sample_split flag, read gcloud bucket to get filenames #
    ############################################################
    # Fixes issue with loading VCF files with wildcard requiring same samples in all VCFs
    if sample_split:
        vcf_files_tmp = []
        for vcf_file in vcf_files:
            if sample_split:
                read_cmd = f"gsutil ls {os.path.join(data_dir, vcf_file)}"
                file_list = subprocess.check_output(shlex.split(read_cmd))

                file_list_clean = str(file_list).replace("\\t", "").split("\\n")
                file_list_clean = [x.replace("b'", "") for x in file_list_clean if x != "'"]

                vcf_files_tmp.extend(file_list_clean)

        vcf_files = vcf_files_tmp

    #############
    # Load VCFs #
    #############
    counter = 1
    row_tmp_fn = os.path.join(out_dir, f"/tmp/{combined_mt_fn.replace('.mt/','')}_rows.ht/")
    for vcf in vcf_files:
        #############################################################################################
        # If MT does not already exist, load in VCF and then write it to disk as a single-sample MT #
        #############################################################################################
        vcf_name = os.path.join(data_dir, vcf)

        vcf_stem = vcf.replace("/","").replace(".vcf", "").replace(".gz", "").replace(".bgz", "").replace("*", "")
        mt_name = os.path.join(out_dir, vcf_stem + f"{test_str}.mt/")

        if (not check_exists(mt_name)) or force:
            logging.info(f'Detected {mt_name} does not exist, importing vcf.')
            # Recode chromosomes/contigs if needed
            if recode is None:
                mt_tmp = hl.import_vcf(
                    vcf_name, force_bgz=force_bgz, call_fields=call_fields,
                    reference_genome=reference_genome, force=force_load)
            else:
                mt_tmp = hl.import_vcf(
                    vcf_name, force_bgz=force_bgz, call_fields=call_fields,
                    reference_genome=reference_genome, contig_recoding=recode, force=force_load)

            # Filter to chr22 only if test flag given
            if test:
                logging.info('Test flag given, filtering to chrom 22.')
                if reference_genome == "GRCh38":
                    chrom_code = "chr22"
                else:
                    chrom_code = "22"

                mt_tmp = mt_tmp.filter_rows(mt_tmp.locus.contig == chrom_code)

            mt_tmp = mt_tmp.checkpoint(mt_name, overwrite=True)
        else:
            logging.info(f"Detected mt of input vcf {vcf} already exists, reading mt directly.")
            mt_tmp = hl.read_matrix_table(mt_name)

        # Annotate matrix table with input file name
        mt_tmp = mt_tmp.annotate_cols(input_file=vcf)
        logging.info('%s imported count: %s' % (vcf, mt_tmp.count()))

        #################################################
        # Combine VCF with main matrix table with union #
        #################################################
        if counter == 1:
            counter += 1
            if save_row_annots:
                row_info = mt_tmp.rows()
                row_info = row_info.checkpoint(row_tmp_fn, overwrite=True)

            mt = mt_tmp

        else:
            counter += 1
            if save_row_annots:
                row_tmp = mt_tmp.rows()
                row_info = row_info.join(row_tmp, how="left")
                row_info = row_info.checkpoint(row_tmp_fn, overwrite=True)

            if chrom_split:
                mt = mt.union_rows(mt_tmp)
            else:
                mt = mt.union_cols(mt_tmp)

        # Annotate matrix table with row info from all the separate VCFs
        if save_row_annots and (len(vcf_files) > 1):
            mt = mt.annotate_rows(**row_info.index(mt.row_key))

        # Checkpoint to avoid running out of memory
        mt = mt.checkpoint(combined_mt_fn, overwrite=True)

    return mt


def check_exists(filename):
    """
    Check if named file exists, adding needed slash if .ht or .mt object
    :param filename:
    :return:
    """
    if filename.startswith("gs://"):
        if filename.endswith(".mt") or filename.endswith(".ht"):
            filename = filename + "/"
        if filename.endswith(".mt/") or filename.endswith(".ht/"):
            filename = filename + "metadata.json.gz" # needs file, not mt/ht dir
        qstat_cmd = ['gsutil', '-q', 'stat', filename]
        exists = subprocess.call(qstat_cmd)

        if exists == 0:
            return True
        else:
            return False

    else:
        return os.path.exists(filename)


def get_path_info(path):  # shamelessly stolen from Pietro
    file_path = os.path.dirname(path)
    basename = os.path.basename(path)
    file_root, file_extension = os.path.splitext(basename)
    return file_path, file_root, file_extension


def return_open_func(f):  # shamelessly stolen from Pietro
    """
    Detects file extension and return proper open_func
    """

    file_path, file_root, file_extension = get_path_info(f)

    if 'bgz' in file_extension:
        open_func = partial(gzip.open, mode='rb')
    elif 'gz' in file_extension:
        open_func = partial(gzip.open, mode='rt')
    else:
        open_func = open
    return open_func


def check_delimiter(f):  # shamelessly stolen from Pietro
    open_func = return_open_func(f)
    with open_func(f) as i:header = i.readline().strip()
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(header)
    if dialect.delimiter != "0":
        return dialect.delimiter
    else:
        return None


def copy_logs_output(log_dir, log_file, plot_dir):
    if not log_dir.endswith("/"):
        log_dir = log_dir + "/"

    datestr = time.strftime("%Y.%m.%d")
    hail_log_name = os.path.join(log_dir, datestr + "_hail_log.txt")
    hl.copy_log(hail_log_name)

    cmd = ['gsutil', 'cp', log_file, log_dir]
    subprocess.call(cmd)
    cmd = ['gsutil', 'cp', '*.html', plot_dir]
    subprocess.call(cmd)
    cmd = ['gsutil', 'cp', '*.pdf', plot_dir]
    subprocess.call(cmd)


def add_secondary(cluster_name, num_secondary, region='europe-west1'):
    """
    Add preemptible nodes on the cluster for applicable (non-shuffle) steps.
    :param cluster_name: Name of the cluster
    :param num_secondary: Number of preemptibles
    :param region: region (default europe-west1)
    :return:
    """
    if cluster_name is not None:
        logging.info(f"Adding {str(num_secondary)} secondary workers to cluster {cluster_name}")

        cmd = shlex.split(f"gcloud dataproc clusters update {cluster_name} --region {region} "
                          f"--num-secondary-workers {str(num_secondary)}")

        subprocess.call(cmd)


def remove_secondary(cluster_name, region='europe-west1'):
    """
    Removes all secondary workers from a cluster.
    :param cluster_name: cluster name
    :param region: region, default 'europe-west1'
    :return:
    """
    if cluster_name is not None:
        logging.info('Removing secondary workers.')

        cmd = shlex.split(f"gcloud dataproc clusters update {cluster_name} --region {region} "
                          f"--num-secondary-workers 0")

        subprocess.call(cmd)


def add_primary(cluster_name, num_primary, region='europe-west1'):
    if cluster_name is not None:
        logging.info(f"Adding {str(num_primary)} primary workers to cluster {cluster_name}")

        cmd = shlex.split(f"gcloud dataproc clusters update {cluster_name} --region {region} "
                          f"--num-workers {str(num_primary)}")

        subprocess.call(cmd)


def remove_primary(cluster_name, region='europe-west1'):
    if cluster_name is not None:
        logging.info('Removing extra primary workers.')

        cmd = shlex.split(f"gcloud dataproc clusters update {cluster_name} --region {region} "
                          f"--num-workers 2")

        subprocess.call(cmd)


def check_counts(mt, args):
    counts = mt.count()
    defined_rows = mt.aggregate_rows(hl.agg.count_where(hl.is_defined(mt.locus)))
    defined_cols = mt.aggregate_cols(hl.agg.count_where(hl.is_defined(mt.s)))
    if (counts[0] != args.start_rows) or (counts[1] != args.start_cols):
        print(f"Error! counts mismatch! Start rows: {args.start_rows}, start cols: {args.start_cols}, count: {counts}")
    else:
        print("counts OK")

    if (defined_rows != args.start_rows) or (defined_cols != args.start_cols):
        print(f"Error! defined mismatch! Start rows: {args.start_rows}, start cols: {args.start_cols}, "
              f"defined rows: {defined_rows}, defined cols: {defined_cols}")
    else:
        print("defined counts OK ")


def tmp_bash(cmd, check=False):
    """
    Creates a temporary bash file to run subprocess command files, useful when you want to have the output piped to
    another file
    :param cmd: command to run
    :param check:
    :return:
    """
    script_file = NamedTemporaryFile(delete=True)

    with open(script_file.name, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(cmd + "\n")
    os.chmod(script_file.name, 0o777)
    script_file.file.close()

    if check:
        subprocess.check_call(script_file.name)
    else:
        subprocess.call(script_file.name, stderr=subprocess.DEVNULL)


def check_vep(mt):
    try:
        test = hl.is_defined(mt.vep)
    except Exception as e:
        logging.error("Error! Input matrix table has not been VEP annotated!")
        logging.error(e)
        exit()


def check_multi_split(mt):
    try:
        test = hl.is_defined(mt.was_split)
    except Exception as e:
        logging.error("Error! Multi-allelic variants must be split before running.")
        logging.error(e)
        exit()


def create_test_dataset(mt, reference_genome, mt_name, out_dir):
    logging.info('Test flag given, filtering to chromosome 22 and chrom X.')

    if reference_genome == "GRCh38":
        chrom_codes = hl.array(["chr22", "chrX"])
    else:
        chrom_codes = hl.array(["22", "X"])

    mt = mt.filter_rows(chrom_codes.contains(mt.locus.contig))

    if mt_name.endswith(".mt/"):
        test_mt = (mt_name.replace(".mt/", "") + "_test.mt").split("/")[-1]
    else:
        test_mt = (mt_name[:-1] + "_test.mt").split("/")[-1]

    mt = mt.checkpoint(os.path.join(out_dir, test_mt), overwrite=True)

    return mt


def get_upper_triangle(list):
    """
    Given a list of items, returns a list where each element is the upper triangle of a pairwise comparison between
    the items in the list.
    :param file_list:
    :return: Returns list of pairs, lenth n choose 2 (where n is length of file list)
    """
    # Get upper triangle of pairwise comparisons
    uppertri = []
    pairs = [[x, y] for i, x in enumerate(list) for j, y in enumerate(list) if i != j]  # upper and lower
    for pair in pairs:  # remove duplicates/lower
        pair.sort()
        if pair not in uppertri:
            uppertri.append(pair)

    return uppertri
