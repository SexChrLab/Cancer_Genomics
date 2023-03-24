#! importing join
from os.path import join

# variant specific read count quantification with GATK ASEreadCounts

# Config

configfile: "somatic_mutation_calling_config.json"

# Directories
EXPRESSION_DIR = "/scratch/eknodel/TESLA/RNAseq/variant_specific_expression/"
vep_path = "~/TESLA_Cancer_Genomics/Cancer_Genomics/02_variant_annotation/variants"
SORTED_BAM_AL_DIR = "/scratch/eknodel/TESLA/RNAseq/sorted_bam/" # path to directory for sorted BAM alignment files

rule all:
	input:
	expand(EXPRESSION_DIR + "{RNA}_gatk_filtered_expression.tsv", EXPRESSION_DIR=EXPRESSION_DIR, RNA=RNA),
        expand(EXPRESSION_DIR + "{RNA}_strelka_unique_rescues.out", EXPRESSION_DIR=EXPRESSION_DIR, RNA=RNA),
        expand(vep_path + "strelka_unique_rescues.vcf.idx", vep_path=vep_path),
        expand(EXPRESSION_DIR + "{RNA}_gatk_filtered_expression.out", EXPRESSION_DIR=EXPRESSION_DIR, RNA=RNA)

rule remove_duplicates:
    input:    
        vars = os.path.join(vep_path, "{sample}_merged_vep_gt.vcf"),
        ref = "/data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    output:
        vars = os.path.join(vep_path, "{sample}_merged_vep_gt_filtered_snps.vcf")
    params: 
        gatk = "/home/eknodel/gatk-4.1.7.0/gatk"
    shell:
        """
        {params.gatk} SelectVariants --V {input.vars} --select-type-to-include SNP -O {output.vars};
        """

rule index_vcf:
    input: 
	vars = os.path.join(vep_path, "{sample}_merged_vep_gt.vcf")
    output:
    	vars = os.path.join(vep_path, "{sample}_merged_vep_gt.vcf.idx")
    params: 
        gatk = "/home/eknodel/gatk-4.1.7.0/gatk"
    shell:
        """
        {params.gatk} IndexFeatureFile -I {input.vars};
        """

rule gatk_asereadcounter:
    input: 
        BAM = lambda wildcards: os.path.join(SORTED_BAM_AL_DIR, config[wildcards.sample]["RNA"] + "_RNA_HISAT2_genome_aligned_sortedbycoord_mkdup_RG.bam"),
    	vars = os.path.join(vep_path, "{sample}_merged_vep_gt.vcf"), 
	vars_idx = os.path.join(vep_path, "{sample}_merged_vep_gt.vcf.idx")
    output:
        vars = os.path.join(EXPRESSION_DIR, "{sample}_filtered_expression.tsv")
    params:
        gatk = "/home/eknodel/gatk-4.1.7.0/gatk",
        ref = "/data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    shell: 
        """
        {params.gatk} ASEReadCounter --output {output.vars} --input {input.BAM} --R {params.ref} --variant {input.vars};  
        """

rule format_output:
    input:
        os.path.join(EXPRESSION_DIR, "{RNA}_filtered_expression.tsv"),
    output:
        os.path.join(EXPRESSION_DIR, "{RNA}_filtered_expression.out"),
    shell:
        """
        cat {input} | awk '{{print $1":"$2-1":"$4":"$5, $6, $7, $8, $12}}' > {output};
        """

