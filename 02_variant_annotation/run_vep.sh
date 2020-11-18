#!/bin/bash
#SBATCH --job-name=vep  # Job name
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=tphung3@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00

PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ perl ~/scratch/Mayo_breast_regional_heterogeneity/01_neoepitope/external_scripts/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl -i /scratch/tphung3/Cancer_Genomics/02_variant_annotation/test.csv --format vcf --cache --assembly GRCh38 --offline --vcf -o /scratch/tphung3/Cancer_Genomics/02_variant_annotation/test_vep.vcf  --force_overwrite --plugin Wildtype --dir_plugins /home/tphung3/.vep/Plugins --symbol --terms SO --plugin Downstream