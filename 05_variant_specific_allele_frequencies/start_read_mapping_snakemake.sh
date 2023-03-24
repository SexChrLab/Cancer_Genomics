#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 1            # number of cores
#SBATCH -t 2-00:00:00   # time in d-hh:mm:ss
#SBATCH --mem=100G
#SBATCH -p general      # partition
#SBATCH -q public       # QOS
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --export=NONE   # Purge the job-submitting shell environment

source activate cancergenomics
PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/per15/site_perl

snakemake --snakefile hisat2_read_mapping.snakefile -j 15 --keep-target-files --rerun-incomplete --cluster "sbatch -N 1 -n 1  -t 2-00:00:00 --mem=100G -p general -q public"

