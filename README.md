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
    4. Mapping
    - Use the Snakefile `map.snakefile`
        1. Prepare reference
        1. Map

3. Variant calling
    - We use three programs for variant calling: VarScan, GATK Mutect2, and Strelka
    - Snakefiles: `varscan.snakefile`, `gatk_mutect2.snakefile`, and `strelka.snakefile`.
    - **Need to edit line 123 of varscan.snakefile to have your own version of VarScan2_format_converter.py. Can be obtained from: https://gist.github.com/PoisonAlien/be1af2a53d5d7bbe2c7a#file-varscan2_format_converter-py** 
    - **Need to edit line 43 of strelka.snakefile to have the appropriate user name on the scratch directory and 92 to appropriate miniconda directory**

## 02_variant_annotation
1. Prepare file as input to the program Variant Effect Predictor (VEP)
- Script `prepare_input_for_vep.py`
    ```
    python prepare_input_for_vep.py --vcf_filenames {/input/full/path/to/vcf/files/} --vep_format_fn {/input/full/path/to/output/vep/file} 
    ```
- This script takes in one or more than one vcf files. If you have more than one VCF file, please separate them by a comma and no space
- GATK file will need to be unzipped using gunzip {filename} before running `prepare_input_for_vep.py`
- If there's just one vcf file as input, the script converts the variant in VCF file format to format that can be used for VEP:
    - Column 1: chromosome name
    - Column 2: Position
    - Column 3 : "."
    - Column 4: reference nucleotide
    - Column 5: alternate nucleotide
    - Column 6, 7, and 8: "."
- If there are more than one vcf files as input, the script selects the variants that are shared in at least 2 vcf files, before converting the variant in VCF file format to the format mentioned above. 
- This script also automatically outputs a file called `number_of_variants_summary.txt` that list the number of variants per vcf file (if there are more than 1 vcf files) and the number of variants that are shared in at least 2 vcf files 

2. Run Variant Effect Predictor (VEP)
- See script `run_vep.sh` for direction on how to run vep. 

## 03_generate_peptide
- We utilize pvacseq (https://anaconda.org/bioconda/pvacseq) to generate peptide. An example command is:
    ```
    pvacseq generate_protein_fasta {/input/to/vep/vcf/} 9 {output}
    ```
- Because the program netMHCpan truncates the name of the transcripts in its output, in order to retain the transcript name information in the netMHCpan outputs, we modify the output from pvacseq using this command:
    ```
    cat {input.peptides} | grep -A 1 ">MT" | sed '/--/d' | sed 's/MT.*ENS//' > {output.peptides_formatted}
    ```

## 04_filtering_peptide
- We are using the following directory structure
- Example: results/Test/HLA-A01:01/9_mers/
- Run netMHCpan example:
    ```
    netMHCpan -a HLA-A01:01 -f ../03_generate_peptides/test_vep_9mers_fmt.txt -BA -s -xls -l 9  -xlsfile results/Test/HLA-A01\:01/9_mers/netmhc.xsl
    ```
- Run netMHCstabpan example:
    ```
    netMHCstabpan -a HLA-A01:01 -f ../03_generate_peptides/test_vep_9mers_fmt.txt -s -xls -l 9  -xlsfile results/Test/HLA-A01\:01/9_mers/netmhcstab.xsl
    ```
  
- Currently (as of November 2020), we are using the thresholds defined for binding affinity, binding stability, and tumor abundance as suggested by Wells et al. (2020) (https://pubmed.ncbi.nlm.nih.gov/33038342/). The repository for the filtering is here: https://github.com/tanyaphung/neoantigens_prioritization. For this example, we do not have RNAseq data, so we are filtering based on binding affinity and binding stability alone. This is how to run the script:
    ```
    python neoepitope_prediction.py --hla_types_fn results/Test/hla.txt --sample_id Test --mers 9 --data_dir results/
    ```
<<<<<<< HEAD
=======
  
## HLA typing
- Can use HLA-LA for HLA typing.
- Refer to the HLA-LA repo (https://github.com/DiltheyLab/HLA-LA) for details on how to get it running.  
>>>>>>> 0bef8179e0749f331d76986fcde2294e5afd23d6
