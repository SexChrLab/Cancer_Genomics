# Context ## TODO:
import os

configfile: "somatic_mutation_calling_config.json"

rule all:
    input:
        expand("fastq_fixed/{sample}_1_fixed.fastq", sample=config["all_samples"]),
        expand("fastq_fixed/{sample}_2_fixed.fastq", sample=config["all_samples"])

rule fix_fastq:
    input:
        fq_1 = lambda wildcards: os.path.join(config["fastq_path"], config[wildcards.sample]["fq_1"]),
        fq_2 = lambda wildcards: os.path.join(config["fastq_path"], config[wildcards.sample]["fq_2"])
    output:
        fq_1 = "fastq_fixed/{sample}_1_fixed.fastq",
        fq_2 = "fastq_fixed/{sample}_2_fixed.fastq"
    shell:
        """
        sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" {input.fq_1} > {output.fq_1};
        sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" {input.fq_2} > {output.fq_2}
        """
