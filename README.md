## Purpose
This library is a collection of tools to analyze and run quality control on exome sequecing data with Hail. The scripts 
are explained in the order I usually run them.

## Import VCF and VEP annotate
`import_vcf_vep_annotate.py` takes one or more VCF files, imports them to Hail matrix table format, combines them, and 
then runs Variant Effect Predictor (VEP) annotation on the files. **Note**: this assumes your VCF files are 
not/no longer split by chromosome, as
it combines the VCF files sample-wise. It's (usually) trivial to combine VCFs split by chromosome with
[bcftools concat](http://www.htslib.org/doc/bcftools.html#concat) first.

To submit the script to a Google Cloud dataproc cluster, use `hailctl dataproc submit` 
([more info](https://hail.is/docs/0.2/cloud/google_cloud.html).) 

### Start your cluster with the `--vep flag!
If you are submitting to a Google Cloud dataproc cluster, for this script you must start your cluster with 
`hailctl dataproc start NAME --vep GENOME_BUILD`

## Exome sequencing data QC
