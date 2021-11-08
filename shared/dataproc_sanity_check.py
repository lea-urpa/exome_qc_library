"""
Functions to check that data inputs to dataproc pipeline are in the same continent as the dataproc cluster.
Author: Lea Urpa, November 2021
"""

import subprocess
import shlex


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
    bucket_describe = f"gsutil ls -L -b {file_url}"

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
