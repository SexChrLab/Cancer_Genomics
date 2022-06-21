# Cancer_Genomics
Assembling pipelines for our cancer genomics work. This includes neoepitope identification, and sex difference analyses.

## Conda environment
- The yml file for the environment (`cancer_genomics_environment.yml`) is located under 00_misc.
- Create the conda environment:
    ```
    conda env create -f cancer_genomics_environment.yml
    ```

## Packages for download:
1. VarScan
    - Download VarScan version 2.3.9 here: https://sourceforge.net/projects/varscan/files/
2. GATK
    - wget https://github.com/broadinstitute/gatk/releases/download/4.1.8.1/gatk-4.1.8.1.zip`
3. Strelka:
    - See: https://github.com/Illumina/strelka
    
## External scripts:
- External scripts are located in the directory `external_scripts`

## 00_misc
- This directory contains miscellaneous files needed to run the program
1. `samples_info.csv`: this is a csv file with 4 columns: sampleID, tissue, type (normal or tumor), and genome type (DNA or RNA).
2. `adapter_sequence.fa`
## 01_somatic_mutation_calling
1. Generate a config file:
    - We generate a config file called `somatic_mutation_calling_config.json` that we will be using for this section. Below is an example of how to run this script. **You need to change the path to the files as appropriate.**
    - Note: This program is written such that normal is always assumed to be listed in the config file before tumor.
        ```
        python generate_config.py --fastq_path /data/CEM/shared/controlled_access/Beauty/
                                  --sample_info /scratch/tphung3/Cancer_Genomics/00_misc/samples_info.csv
                                  --ref_dir /data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome
                                  --ref_basename GRCh38_full_analysis_set_plus_decoy_hla
                                  --varscan_path /home/tphung3/softwares/VarScan.v2.3.9.jar
                                  --gatk_path /home/tphung3/softwares/gatk-4.1.8.1/gatk 
                                  --strelka /home/tphung3/softwares/miniconda3/envs/cancer/share/strelka-2.9.10-0/bin/configureStrelkaSomaticWorkflow.py
                                  --bam_readcount /home/tphung3/softwares/bam-readcount
                                  --perl_fp_filter /home/tphung3/softwares/fpfilter-2.pl
        ```
2. Quality control
    - Use the Snakefile `quality_control.snakefile`:
        1. Raw QC
        1. Trim
        1. Trimmed QC
    - **You need to edit lines 5-6 of the snakefile to point to the appropriate path to the adaptor sequence and perl lib for your system.**

3. Mapping
    - Use the Snakefile `map.snakefile`
        1. Prepare reference
        2. Map
        3. Mark duplicates
        4. Index

4. Variant calling
    - We use three programs for variant calling: VarScan, GATK Mutect2, and Strelka
    - Snakefiles: `varscan.snakefile`, `gatk_mutect2.snakefile`, and `strelka.snakefile`.
    - **Need to edit line 43 of strelka.snakefile to have the appropriate user name on the scratch directory and 92 to appropriate miniconda directory**

## 02_variant_annotation

Snakefile `VEP_PVACseq.snakefile` performs the following steps:  

1. Prepare file as input to the program Variant Effect Predictor (VEP)
    - GATK file is unzipped
    - Runs the script `prepare_input_for_vep.py`
    - As currently written, the snakemake automatically uses the outputs of GATK and Strelka only and takes the overlap of these two programs. 
    - Other available options are:  
        - One vcf file can be provided as input, at which pointthe script converts the variant in VCF file format to format that can be used for VEP:
            - Column 1: chromosome name
            - Column 2: Position
            - Column 3 : "."
            - Column 4: reference nucleotide
            - Column 5: alternate nucleotide
            - Column 6, 7, and 8: "."
        - Additional VCF files can be merged by adding them to the list with a comma and no space. Note that the program will select the overlap of any two of the VCF files provided.  
    - This script also automatically outputs a file called `number_of_variants_summary.txt` that list the number of variants per vcf file (if there are more than 1 vcf files) and the number of variants that are shared in at least 2 vcf files 

2. Runs Variant Effect Predictor (VEP)
    - Returns annotated variants 
    - Uses downstream pluggin to eliminate variants that are downstream of a frameshift mutation
    - Uses wildtype pluggin to return the corresponding wildtype peptide as well

3. Runs PVAC-Seq (https://anaconda.org/bioconda/pvacseq) to generate peptide. 
    - The default command is:
        ```
        pvacseq generate_protein_fasta {/input/to/vep/vcf/} 17 {output}
        ```
    - This will generate 17mer peptides but changing the 17 can change the default for the generated peptides. Generating 17mers is recommended as this allows the following steps to create 9mer peptides in which the mutation of interest is in every possible position

## 03_calculate_neoantigen_binding
We assess neoantigen binding in two ways. First, we use netMHCpan to calculate the dissociation constant and then we use the netMHCstabpan to calculate the binding stability

**Note, currently the HLA types are assumed to be present in the config file, but the config file created here does not automatically add them, this will be fixed soon**

All of the following steps are performed with `netMHCpan-snakemake.py`:

1. Formatting input

2. netMHCpan for calculation of dissociation constant

3. netMHCstabpan for calculation of stability

## 04_expression
- We utilize salmon for expression quantification (https://salmon.readthedocs.io/en/latest/salmon.html)
- Need to have a salmon index and then change lines 13, 16, and 17. **To do: Add these directories to the config file to reduce hard coding**

## 05_filtering_peptide

**To be updated**

- Currently (as of November 2020), we are using the thresholds defined for binding affinity, binding stability, and tumor abundance as suggested by Wells et al. (2020) (https://pubmed.ncbi.nlm.nih.gov/33038342/). The repository for the filtering is here: https://github.com/tanyaphung/neoantigens_prioritization. For this example, we do not have RNAseq data, so we are filtering based on binding affinity and binding stability alone. This is how to run the script:
    ```
    python neoepitope_prediction.py --hla_types_fn results/Test/hla.txt --sample_id Test --mers 9 --data_dir results/
    ```
  
## HLA typing

**To be updated**
- Can use HLA-LA for HLA typing.
- Refer to the HLA-LA repo (https://github.com/DiltheyLab/HLA-LA) for details on how to get it running.  
