import os

configfile: "somatic_mutation_calling_config.json"

rule all:
    input:
        expand("strelka/{subject}/runWorkflow.py", subject=config["all_subjects"]),
        expand("strelka/{subject}/results/variants/somatic.snvs.vcf.gz", subject=config["all_subjects"]),
        expand("strelka/{subject}/results/variants/somatic.indels.vcf.gz", subject=config["all_subjects"])

rule config:
    input:
        normal_bam = lambda wildcards: os.path.join(
			"processed_bams/", config[wildcards.subject]["normal"] + "." + config["ref_basename"] + ".sorted.chr21.bam"),
        tumor_bam = lambda wildcards: os.path.join(
			"processed_bams/", config[wildcards.subject]["tumor"] + "." + config["ref_basename"] + ".sorted.chr21.bam"),
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa")
    output:
        "strelka/{subject}/runWorkflow.py"
    params:
        strelka = config["strelka"],
        run_dir = "strelka/{subject}"
    shell:
        """
        {params.strelka} --normalBam {input.normal_bam} --tumorBam {input.tumor_bam} --referenceFasta {input.ref} --runDir {params.run_dir}
        """

rule run:
    input:
        normal_bam = lambda wildcards: os.path.join(
			"processed_bams/", config[wildcards.subject]["normal"] + "." + config["ref_basename"] + ".sorted.chr21.bam"),
        tumor_bam = lambda wildcards: os.path.join(
			"processed_bams/", config[wildcards.subject]["tumor"] + "." + config["ref_basename"] + ".sorted.chr21.bam"),
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa")
    output:
        snvs = "strelka/{subject}/results/variants/somatic.snvs.vcf.gz",
        indels = "strelka/{subject}/results/variants/somatic.indels.vcf.gz"
    params:
        run = "/scratch/tphung3/Cancer_Genomics/01_somatic_mutation_calling/strelka/{subject}/runWorkflow.py"
    shell:
        """
        {params.run} -m local -j 20
        """
