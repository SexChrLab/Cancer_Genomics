# Setting up filesnames here:
from os.path import join
import os

configfile: "somatic_mutation_calling_config.json"

rule all:
    input:
        expand("variants/{sample}_merged.vcf", sample=config["all_subjects"]), # combine all mutations
        expand("variants/{sample}_merged_vep.vcf", sample=config["all_subjects"]), # run VEP
        expand("peptides/{sample}_vep.17.peptide", sample=config["all_subjects"]) # run pvacseq

rule gunzip:
    input:
        gatk = os.path.join("../01_somatic_mutation_calling/gatk_mutect2/{sample}.somatic.filtered.pass.vcf.gz")
    output:
        gatk = os.path.join("../01_somatic_mutation_calling/gatk_mutect2/{sample}.somatic.filtered.pass.vcf")
    shell: 
        """
        gunzip {input.gatk}
        """

rule generate_input:
    input:
        gatk = os.path.join("../01_somatic_mutation_calling/gatk_mutect2/{sample}.somatic.filtered.pass.vcf"),
        strelka = os.path.join("../01_somatic_mutation_calling/strelka/{sample}/results/variants/somatic.snvs.pass.vcf"),
        strelka_indels = os.path.join("../01_somatic_mutation_calling/strelka/{sample}/results/variants/somatic.indels.pass.vcf") 
    output:
        vep = os.path.join("variants/{sample}_merged.vcf")
    shell: 
        """
        python prepare_input_for_vep.py --vcf_filenames \
        {input.gatk},{input.strelka},{input.strelka_indels} \
        --vep_format_fn {output.vep}
        """

rule sort:
    input:
        vcf = os.path.join("variants/{sample}_merged.vcf")
    output:
        vcf = os.path.join("variants/{sample}_merged_sorted.vcf")
    shell:
        """
        sort -V {input.vcf} > {output.vcf}
        """

rule run_vep:
    input:
        os.path.join("variants/{sample}_merged_sorted.vcf")
    output:
        os.path.join("variants/{sample}_merged_vep.vcf")
    shell:
        """
        vep -i {input} --format vcf --assembly GRCh38 --cache --dir_cache ~/external_scripts --offline --vcf -o {output} --force_overwrite --plugin Wildtype --symbol --terms SO --plugin Downstream
        """ 

rule generate_fasta:
    input:
        os.path.join("variants/{sample}_merged_vep.vcf")
    output:
        len_21 = os.path.join("peptides/{sample}_vep.17.peptide"),
    shell:
        """
        pvacseq generate_protein_fasta {input} 17 {output.len_21};
        """
