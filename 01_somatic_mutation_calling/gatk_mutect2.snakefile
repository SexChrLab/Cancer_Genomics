import os

configfile: "somatic_mutation_calling_config.json"

rule all:
    input:
        expand(os.path.join("gatk_mutect2/", "{subject}.somatic.vcf.gz"), subject=config["all_subjects"]),
        expand(os.path.join("gatk_mutect2/", "{subject}.somatic.filtered.vcf.gz"), subject=config["all_subjects"])

rule tumor_with_matched_normal:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        normal_bam = lambda wildcards: os.path.join(
			"processed_bams/", config[wildcards.subject]["normal"] + "." + config["ref_basename"] + ".sorted.chr21.bam"),
        tumor_bam = lambda wildcards: os.path.join(
			"processed_bams/", config[wildcards.subject]["tumor"] + "." + config["ref_basename"] + ".sorted.chr21.bam")
    output:
        os.path.join("gatk_mutect2/", "{subject}.somatic.vcf.gz")
    params:
        gatk = config["gatk_path"],
        sm = lambda wildcards: config[wildcards.subject]["normal"]
    shell:
        """
        {params.gatk} Mutect2 -R {input.ref} -I {input.tumor_bam} -I {input.normal_bam} -normal {params.sm} -O {output}
        """

rule filter:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        unfiltered = os.path.join("gatk_mutect2/", "{subject}.somatic.vcf.gz")
    output:
        filtered = os.path.join("gatk_mutect2/", "{subject}.somatic.filtered.vcf.gz")
    params:
        gatk = config["gatk_path"]
    shell:
        """
        {params.gatk} FilterMutectCalls -R {input.ref} -V {input.unfiltered} -O {output.filtered}
        """
