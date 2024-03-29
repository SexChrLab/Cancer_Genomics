# Cancer_Genomics

**Last updated:** 6/7/2023

**Updated by:** Elizabeth Borden

**Contact information:** knodele@email.arizona.edu

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
    - If sex specific read mapping is desired, generate the config separately for read mapping in the folder sex_specific_read_mapping with the following command
        ```
        python generate_config.py --fastq_path /data/CEM/shared/controlled_access/Beauty/
                                  --sample_info /scratch/tphung3/Cancer_Genomics/00_misc/samples_info.csv
                                  --female_ref_dir /data/CEM/shared/public_data/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC/
                                  --female_ref_basename GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded
                                  --male_ref_dir /data/CEM/shared/public_data/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC/
                                  --male_ref_basename GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded
                                  --ref_basename GCA_009914755.4_CHM13_T2T_v2.0
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
    - For sex specific mapping, use the Snakefile `sex_specific_read_mapping/map.snakefile`, this will perform the same steps as outlines above and should seemlessly integrate back into the main pipeline

4. Variant calling
    - We use three programs for variant calling: VarScan, GATK Mutect2, and Strelka
    - Snakefiles: `varscan.snakefile`, `gatk_mutect2.snakefile`, and `strelka.snakefile`.
    - **Need to edit line 43 of strelka.snakefile to have the appropriate user name on the scratch directory and 92 to appropriate miniconda directory**
    - **TO DO: check if we need to add a sex specific option for the variant calling**

## 02_variant_annotation

Snakefile `VEP_PVACseq.snakefile` performs the following steps:  

1. Prepare file as input to the program Variant Effect Predictor (VEP)
    - GATK file is unzipped
    - Runs the script `merge_Mutect2_Strelka2_variants.py`
    - As currently written, the snakemake automatically uses the outputs of GATK and Strelka only and takes the overlap of these two programs. **Update 6/7:** A bug was discovered in the original script used here and we determined that multinucleotide variants were being discarded. We have therefore created a new script that rescues these multinuecleotide variants. However, as of 6/7 this only works for the overlap of GATK Mutect2 and Strelka2. The original file `prepare_input_for_vep.py` is still provided, but should not be used for overlapping variant callers
    - File overlap is also formatted in a VCF file format that can be used for VEP:
            - Column 1: chromosome name
            - Column 2: Position
            - Column 3 : "."
            - Column 4: reference nucleotide
            - Column 5: alternate nucleotide
            - Column 6, 7, and 8: "."
    - This script also automatically outputs a file called `number_of_variants_summary.txt` that list the number of variants per vcf file and the number of variants that are shared in at least 2 vcf files 

2. Runs Variant Effect Predictor (VEP) - **Note:** VEP can be difficult to download. Detailed instructions included below. 
    - Returns annotated variants 
    - Uses downstream pluggin to eliminate variants that are downstream of a frameshift mutation
    - Uses wildtype pluggin to return the corresponding wildtype peptide as well
    - **NOTE: Currently not compatable with the CHM13 reference genome, working on this**

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

## 05_variant_specific_allele_frequency

1. Read mapping of RNAseq data peformed with Hisat2

2. Variant specific allele frequency determined with GATK ASEReadCounts
  
## HLA typing

**To be updated**
- Can use HLA-LA for HLA typing.
- Refer to the HLA-LA repo (https://github.com/DiltheyLab/HLA-LA) for details on how to get it running.  

## Instructions for installing the Variant Effects Predictor

### With Conda

1. Create VEP environment
```
conda create --name vep_env
```

2. Install VEP
```
conda install -c bioconda ensembl-vep -n vep_env
```

3. Download cache
```
cd $HOME/.vep
curl -O https://ftp.ensembl.org/pub/release-109/variation/indexed_vep_cache/homo_sapiens_vep_109_GRCh38.tar.gz
tar xzf homo_sapiens_vep_109_GRCh38.tar.gz
```

4. Add plugins to path - move from the external scripts directory of this repository to ~/.vep/Plugins

### Directly from GitHub Source

* See Ensembl directions: https://uswest.ensembl.org/info/docs/tools/vep/script/vep_download.html

1. Before installing VEP, download or load Perl 5 (eg. perl/5.26.1 module in ASU Agave)
2. Navigate to local directory where you store software
3. Clone ensembl-vep from Github and enter the directory
    ```
    git clone https://github.com/Ensembl/ensembl-vep.git
    cd ensembl-vep
    ```
4. Run installation script which will get the required Perl libraries needed to run VEP
    ```
    perl INSTALL.pl
    ```
5. Select to cache annotation information for any genome references you align to regularly and select to install the Wildtype and Downstream plugins plus any others that sound interesting 
6. Export your install directory to your perl 5 path: 
    ```
    export PERL5LIB=${PERL5LIB}:/your_install_path/ensembl-vep/
    ```
