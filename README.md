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

## Steps in QC process
The exome QC pipeline runs through a number of steps for weeding out failing samples and variants
from your dataset. In general, however, *no samples or variants are excluded*- they are only marked with annotations
that they are failing, and only excluded before particular steps that require them to be removed.

A named checkpoint is written after each of these steps in the output directory specified, and you can restart the 
pipeline after a given checkpoint with the flag `--checkpoint`.

General options for the pipeline:
```
-mt                         Matrix table to run QC on, output from import_vcf_vep_annotate.py
--checkpoint                Load pipeline from specified checkpoint
--reference_genome          Of either GRCh37 or GRCh38 (case sensitive)
--test                      Restricts to chromosome 22 + chromosome X for testing the pipeline
--overwrite_checkpoints     Overwrites previous checkpoints, default true
--run_king                  Stops the pipeline after checkpoint 4, to run King outside of Hail
--num_preemptible_workers   Number of preemptible workers for scaling the cluster in steps which it is allowed
--cluster_name              Name of the Google cloud dataproc cluster the script is submitted to, for scaling
--region                    Region name for scaling in the pipeline
```
Input and output directories for the pipeline:
```
--scripts_dir               Location of this library in the cloud (e.g. gs://my_bucket/exome_qc_library)
--log_dir                   Directory to write logs to (e.g. gs://my_bucket/logs)
--out_dir                   Directory to write output files to (e.g. gs://my_bucket/qc_output)
--out_name                  Output name
```

### Step 1: Annotate Samples
Step one takes a list of files given with `--samples_annotation_files` (comma separated, if more than one) and annotates
the columns of the matrix table with all the columns in the annotation file.

Options for samples annotation:
```
--samples_annotation_files  File names (e.g. gs://my_bucket/phenos/samples_info.txt), comma separated if more than one
--pheno_col                 Column name with case/control status in true/false format.
--samples_col               Column name of sample IDs in sample annotation files (same in all files!)
--samples_delim             Delimiter for sample annotation files, default tab (same in all files!)
--samples_miss              Missing character in sample annotation files, default NA (same in all files!)
--bam_metadata              File containing bam metadata (chimera + contamination info). If not given, assumed the data
                            are in one of the sample annotation files.
--bam_sample_col            Column name of sample IDs in bam metadata file.
--bam_delim                 Delimiter for bam metadata file, default tab.
--bam_miss                  Missing character in bam metadata file, default NA.
--fam_id                    Column name for family ID, if given in sample annotation files
--pat_id                    Column name for paternal ID, if given in sample annotation files
--mat_id                    Column name for maternal ID, if given in sample annotation files
```

### Step 2: Samples removal
Samples can be removed at this step from either a given text file with one sample name per line, or substring(s) in
the sample name to search for.

Sample removal options:
```
--sample_removal_list       List of samples to remove, one sample name per line (e.g. gs://my_bucket/remove_samples.txt)
--sample_removal_strings    substrings, comma separated if more than 1. Samples will be removed if their ID contains 
                            this string.
```

### Step 3: Low-pass variant QC
Variants and genotypes are first removed with lower thresholds, prior to samples QC, to get rid of the 'worst' variants 
before looking for samples that fail QC.

Thresholds specific to low pass variant QC:
```
--low_pass_p_hwe            Hardy-Weinberg p value threshold, default 1e-9
--low_pass_min_call_rate    Minimum variant call rate (number of individuals with a valid call at that variant), 
                            default 0.8
```
Thresholds that are the same for both low-pass and final genotype QC:
```
--snp_qd                    Minimum quality by depth for a snp, default 2
--indel_qd                  Minimum quality by depth for an indel, default 3
--ab_allowed_dev_het        Percentage of het GT calls that must be in allelic balance, default 0.8
``` 

Genotype thresholds (same for both low-pass and final genotype QC):
``` 
--min_dp                    Minimum read depth for a genotype, default 10
--min_gq                    Minimum GQ for a genotype (actually min GQ for hom ref reads, min PL[ref] for het and hom 
                            alt genotypes), default 20
--min_het_ref_reads         Minimum percentage of reference reads for a het genotype, default 0.2
--max_het_ref_reads         Maximum percentage of reference reads for a het genotype, default 0.8
--min_hom_ref_ref_reads     Minimum percentage of reference reads for a hom ref genotype, default 0.9
--max_hom_alt_ref_reads     Maximum percentage of reference reads for a hom alt genotype, default 0.1

```
### Step 4: LD prune and MAF filter dataset for relatedness calculation with King
This step is a step of its own because it takes a while, so it's useful to have a checkpoint after it. This step 
removes variants and genotypes failing low pass variants QC, MAF filters the dataset, LD prunes it, downsamples it
if the variant count is still more that 80,000 variants.

Options for LD pruning and MAF filtering for relatedness:
```
--r2                        R^2 correlation cutoff for LD pruning, default 0.2
--bp_window_size            Window size to look for LD pruning correlations, in base pairs. Default 500,000 bp.
--ind_maf                   MAF cutoff for excluding rare variants, default 0.001
```

### Step 5: Find related individuals
This step either exports the MAF filtered, LD pruned, and downsampled dataset to Plink format to run King
relatedness calculations, (see get_relatedness.py), or uploads a list of related individuals to remove for further
steps where it's necessary to use only unrelated individuals.

Options for finding related individuals:
```
--run_king                  If this flag is given, exports file to Plink and exits pipeline.
--relatives_removal_file    List of samples to remove to give maximum unrelated set (output of get_relatedness.py).
                            If run_king is not given, this list must be given.
```

### Step 6: Find population outliers
After annotating related individuals from a file, those individuals (and failing samples and genotypes) are removed
from the dataset to find population outliers. Principal components are calculated on the data, and any samples that are
more than 4 standard deviations outside the mean on PC1 or PC2 are marked as population outliers and removed from the
dataset. This is repeated until no samples are more that 4 standard deviations outside the mean on PC1 or PC2. 
Population outlier samples are not removed from the data, just marked as population outliers.

Options for population outliers:
``` 
--pop_sd_threshold          Number of standard deviations outside the mean to call a pop outlier, default 4.
--pca_plot_annotations      Sample column names to annotate the PCA plots made at each round, commma separated if more
                            than one.
```

### Step 7: Impute sex
Impute the sex of the samples and annotate samples and variants with sex-aware call rates for chromosome X variants.

Options for imputing sex:
``` 
--female_threshold          F statistic threshold for calling a female, default 0.4
--male_threshold            F statistic threshold for calling a male, default 0.8
```

### Step 8: Samples QC
Finds samples failing on global measures (chimera percentage, contamination percentage) or by deviation from the mean
in their cohort (TiTv ratio, insertion/deletion ratio, het/hom var ratio, number of singletons).

Options for samples QC:
``` 
--chimeras_max              Max fraction of chimeric reads allowed for a sample, default 0.05
--chimeras_col              Column in sample annotations giving chimeras info
--contamination_max         Max fraction of contaminated reads allowed for a sample, default 0.05
--contamination_col         Column in sample annotations giving contamination info
--batch_col_name            If given, separates data into batches to calculate TiTV, insertion/deletion, het/homvar
                            ratios, and n singletons
--sampleqc_sd_threshold     Number of standard devations from the mean a value can be for TiTv, indel, het/homvar,
                            n singleton values, default 4
```

### Step 9: Final variant QC
Performs final variant QC, with same thresholds as in step 3, except the following:

``` 
--final_p_hwe               Hardy-Weinberg p value, default 1e-6
--final_min_call_rate       Minimum call rate for a variant, default 0.9
```

### Step 10: Filter variants by phenotype
If you have case and control data coming from two different sequencing centers and/or exome capture kits, it could be
that the call rate or allelic balance filters may be passed over all for a variant, but not when you look specifically
at cases or controls. That's key- if a variant is poorly called only in controls, you might get a false positive
when doing analysis. So if `--pheno_col` is given, call rates and allelic balances are calculated separately for cases
and controls, and if either is failing the variant is marked as failing QC.

Options for filtering variants by phenotype:
``` 
--pheno_col                 Column giving case-control status as true/false
--pheno_call_rate           Minimum variant call rate, calculated separately for cases and controls, default 0.95
```
The threshold for fraction of het genotypes in allelic balance for a variant is given by the variant QC threshold,
`--ab_allowed_dev_het`.

### Step 11: Calculate final PCs
This step removes failing variants, genotypes, and samples, and samples that are population outliers. It then also 
removes samples that are relatives, and calculates final principal components to use as covariates in further analyses.
After calculating the PCS on unrelated individuals, the PC loadings are used to project principal components for the 
relatives, so all samples that are passing QC and not population outliers have principal components annotated.

Options for PC calculation:
``` 
--pc_num                    Number of principal components to calculate on QCd datset, default 10
```

## Example Submission Scripts
### Part One- up to King relatedness calculation
King relatedness calculation is not supported in Hail, so we first run the first part of the pipeline and export the
data in Plink format to run King.

``` 
# Start Hail cluster
hailctl dataproc start mycluster --max-idle 10m --scopes cloud-platform

# Submit script
hailctl dataproc submit mycluster /path/to/exome_qc_library/exome_qc.py \
-mt gs://my_bucket/mydata_vep_annotated_mt.mt/ \
--run_king --cluster_name mycluster \
--out_name mydatat_qcd \
--out_dir gs://my_bucket/qc_output/ \
--log_dir gs://my_bucket/qc_output/logs \
--scripts_dir gs://my_bucket/exome_qc_library \
--samples_annotation_files gs://my_bucket/phenotypes/sample_info.txt \
--chimeras_col X.CHIMERA --contamination_col X.CONTAMINATION \
--samples_col SAMPLE_ID_IN_VCF --samples_miss NA  \
--batch_col_name batch_cohort --pheno_col is_case \
--fam_id fam_ID --pat_id pat_ID --mat_id mat_ID \
--sample_removal_strings COHORT1 \
--sample_removal_list gs://my_bucket/phenotypes/bad_samples.txt \
--pca_plot_annotations COHORT,CAPTURE 
```

