#!/bin/bash
#SBATCH --job-name=vep  # Job name
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=tphung3@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00

source activate var_call_env

PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/per15/site_perl/ vep -i /scratch/eknodel/Cancer_Genomics/02_variant_annotation/Combined_output_patient12.vep.in --format vcf --database --assembly GRCh38 --vcf -o /scratch/eknodel/Cancer_Genomics/02_variant_annotation/Combined_output_patient12.vep --force_overwrite  --symbol --plugin Wildtype --terms SO --plugin Downstream
