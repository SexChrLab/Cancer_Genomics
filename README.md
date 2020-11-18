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
1. `samples_info.csv`: this is a csv file with 3 columns: sampleID, tissue, and type (normal or tumor). TODO: add DNA or RNA information.
2. `adapter_sequence.fa`
## 01_somatic_mutation_calling
1. Generate a config file:
    - We generate a config file called `somatic_mutation_calling_config.json` that we will be using for this section
        ```
        python generate_config.py
        ```
    - Before running: edit the information in the python script to reflect the accurate information for your samples
2. Quality control
    - Use the Snakefile `quality_control.snakefile`:
        1. Raw QC
        1. Trim
        1. Trimmed QC
    4. Mapping
    - Use the Snakefile `map.snakefile`
        1. Prepare reference
        1. Map

3. Variant calling
    - We use three programs for variant calling: VarScan, GATK Mutect2, and Strelka
    - Snakefiles: `varscan.snakefile`, `gatk_mutect2.snakefile`, and `strelka.snakefile`.

## 02_variant_annotation
1. Prepare file as input to the program Variant Effect Predictor (VEP)