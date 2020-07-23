import os

configfile: "somatic_mutation_calling_config.json"

rule all:
    input:
        expand("pileups/{sample}.pileup", sample=config["all_samples"]), #run bam_pileup
        expand("intermediate_files/{subject}.varscan.snp", subject=config["all_subjects"]), #run VarScan
        expand("intermediate_files/{subject}.varscan.indel", subject=config["all_subjects"]), #run VarScan

rule bam_pileup: #for both normal and tumor
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        bam = os.path.join("processed_bams/{sample}." + config["ref_basename"] + ".sorted.bam")
    output:
        pileup = "pileups/{sample}.pileup"
    threads: 4
    shell:
        """
        samtools mpileup -f {input.ref} {input.bam} > {output.pileup}
        """

rule run_varscan:
    input:
        normal_pileup = lambda wildcards: os.path.join(
			"pileups/", config[wildcards.subject]["normal"] +  ".pileup"),
        tumor_pileup = lambda wildcards: os.path.join(
			"pileups/", config[wildcards.subject]["tumor"] + ".pileup")
    output:
        snp = "intermediate_files/{subject}.varscan.snp",
        indel = "intermediate_files/{subject}.varscan.indel"
    params:
        varscan = config["varscan_path"],
        basename = "intermediate_files/{subject}.varscan"
    threads: 4
    shell:
        "java -jar {params.varscan} somatic {input.normal_pileup} {input.tumor_pileup} {params.basename} –min-coverage 10 –min-var-freq 0.08 –somatic-p-value 0.05"
