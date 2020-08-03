## Purpose
This library is a collection of tools to analyze and run quality control on exome sequecing data with Hail.

## Script order
Generally, the order you'd want to run the scripts is first importing the data from VCF format, then annotate with VEP 
and save the data as a matrix table (as this step takes a long time) with *import_vcf_vep_annotate.py*. 
