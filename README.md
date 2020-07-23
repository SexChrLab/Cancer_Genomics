# Cancer_Genomics
Assembling pipelines for our cancer genomics work. This includes neoepitope identification, and sex difference analyses.

## Conda environment

## Packages for download:
1. VarScan
  1. Download VarScan version 2.3.9 here: https://sourceforge.net/projects/varscan/files/
2. GATK
  1. `wget https://github.com/broadinstitute/gatk/releases/download/4.1.8.1/gatk-4.1.8.1.zip`
## 00_misc
- This directory contains miscellaneous files needed to run the program
1. `samples_info.csv`: this is a csv file with 3 columns: sampleID, tissue, and type (normal or tumor). #TODO: add DNA or RNA information.
2. `adapter_sequence.fa`
## 01_somatic_mutation_calling
1. Generate a config file:
    1. We generate a config file called `somatic_mutation_calling_config.json` that we will be using for this section
    1. `python generate_config.py`
    1. TODO before running: edit the information in the python script to reflect the accurate information for your samples
2. Edit the fastq files
- If the fastq files are downloaded from SRA, or if you encounter this message when running bwa: `paired reads have different names`
- Use the Snakefile `fix_fastq_from_sra.snakefile`
3. Quality control
- Use the Snakefile `quality_control.snakefile`:
    1. Raw QC
    1. Trim
    1. Trimmed QC
4. Mapping
- Use the Snakefile `map.snakefile`
    1. Prepare reference
    1. Map
