#!/bin/bash
#SBATCH --job-name=Gatk  # Job name
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=eknodel@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 168:00:00
#SBATCH --mem=40000
#SBATCH -q wildfire
#SBATCH -p wildfire

newgrp combinedlab

source activate var_call_env

module load picard/2.18.3
module load bcftools/1.9
PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/per15/site_perl

snakemake --snakefile VEP_PVACseq.snakefile -j 71 --keep-target-files --rerun-incomplete --cluster "sbatch -q wildfire -p wildfire -n 1 -c 1 --mem=50000 -t 96:00:00"
