# Purpose
This library is a collection of tools to analyze and run quality control on exome sequecing data with Hail. The scripts 
are explained in the order I usually run them.

# Prerequisites
[Install Hail on your computer](https://hail.is/docs/0.2/getting_started.html),
 even if you are running in the cloud. This will also install the `hailctl` library 
as well.

Clone this library to your laptop. If you are running in the cloud with `dataproc`, you must also `gsutil cp` this 
library to the cloud as well.

# Import VCF and VEP annotate
`import_vcf_vep_annotate.py` takes one or more VCF files, imports them to Hail matrix table format, combines them, and 
then runs Variant Effect Predictor (VEP) annotation on the files. **Note**: this assumes your VCF files are 
not/no longer split by chromosome, as
it combines the VCF files sample-wise. It's (usually) trivial to combine VCFs split by chromosome with
[bcftools concat](http://www.htslib.org/doc/bcftools.html#concat) first.

To submit the script to a Google Cloud dataproc cluster, use `hailctl dataproc submit` 
([more info](https://hail.is/docs/0.2/cloud/google_cloud.html).) 

## Start your cluster with the `--vep flag!
If you are submitting to a Google Cloud dataproc cluster, for this script you must start your cluster with 
`hailctl dataproc start NAME --vep GENOME_BUILD`

## Example submission
In this example submission, the `import_vcf_vep_annotate.py` script is submitted via `hailctl dataproc submit` from a
local copy of this repo to the cloud, pointing to where the copy of the repo exists
in the bucket for the script to find the helper functions.

```
# Start Hail cluster
hailctl dataproc start vep-test --vep GRCh37

# Submit script
hailctl dataproc submit vep-test \
/local/path/exome_qc_library/import_vcf_vep_annotate.py \
--vcf sampleset_1.vcf.gz,sampleset_2.vcf.gz \
--out_file combined_samples_vep_annotated \
--scripts_dir gs://my-bucket/exome_qc_library \
--log_dir gs://my-bucket/logs/ \
--data_dir gs://my_bucket/input_vcfs/ \
--force_bgz --test
```

# Exome sequencing data QC

## 

## Example Submission Scripts
### Part One- up to King relatedness calculation
King relatedness calculation is not supported in Hail,