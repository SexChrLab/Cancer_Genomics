#!/bin/bash
#SBATCH --job-name=variant_expression  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p wildfire
#SBATCH -q wildfire
#SBATCH -t 96:00:00
#SBATCH --mem=1024

newgrp combinedlab

source activate cancergenomics
module load bcftools/1.9
PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/per15/site_perl

snakemake --snakefile variant_expression_quantification.snakefile -j 15 --keep-target-files --rerun-incomplete --cluster "sbatch -n 1 -c 8 -t 96:00:00 -p wildfire -q wildfire"
