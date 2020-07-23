import os

configfile: "somatic_mutation_calling_config.json"

rule all:
    input:
        expand("pileups/{sample}.pileup", sample=config["all_samples"]), #run bam_pileup
        expand("intermediate_files/{subject}.varscan.snp", subject=config["all_subjects"]), #run VarScan
        expand("intermediate_files/{subject}.varscan.indel", subject=config["all_subjects"]), #run VarScan
        expand("varscan/{subject}.varscan.variants.filter.pass") #filter

rule bam_pileup: #for both normal and tumor
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        bam = os.path.join("processed_bams/{sample}." + config["ref_basename"] + ".sorted.chr21.bam")
    output:
        pileup = "varscan/pileups/{sample}.pileup"
    threads: 4
    shell:
        """
        samtools mpileup -f {input.ref} {input.bam} > {output.pileup}
        """

rule run_varscan:
    input:
        normal_pileup = lambda wildcards: os.path.join(
			"varscan/pileups/", config[wildcards.subject]["normal"] +  ".pileup"),
        tumor_pileup = lambda wildcards: os.path.join(
			"varscan/pileups/", config[wildcards.subject]["tumor"] + ".pileup")
    output:
        snp = "varscan/{subject}.varscan.snp",
        indel = "varscan/{subject}.varscan.indel"
    params:
        varscan = config["varscan_path"],
        basename = "intermediate_files/{subject}.varscan"
    threads: 4
    shell:
        "java -jar {params.varscan} somatic {input.normal_pileup} {input.tumor_pileup} {params.basename} –min-coverage 10 –min-var-freq 0.08 –somatic-p-value 0.05"

rule isolate_calls_by_type_and_confidence:
    input:
        varscan_snp = "varscan/{subject}.varscan.snp"
    output:
        varscan_snp = "varscan/{subject}.varscan.snp.somatic.hc"
    params:
        varscan = config["varscan_path"]
    shell:
        """
        java -jar {params.varscan} processSomatic {input.varscan_snp}
        """

rule somatic_filter:
    input:
        snp_somatic_hc = "varscan/{subject}.varscan.snp.somatic.hc",
        indel = "varscan/{subject}.varscan.indel"
    output:
        snp_somatic_hc_filter = "varscan/{subject}.varscan.snp.somatic.hc.filter"
    params:
        varscan = config["varscan_path"]
    shell:
        """
        java -jar {params.varscan} somaticFilter {input.snp_somatic_hc} -indel-file {input.indel} -output-file {output.snp_somatic_hc_filter}
        """

rule convert_to_bed_fmt:
    input:
        snp_somatic_hc_filter = "varscan/{subject}.varscan.snp.somatic.hc.filter"
    output:
        snp_somatic_hc_filter_bed = "varscan/{subject}.varscan.snp.somatic.hc.filter.bed"
    shell:
        """
        awk -F "\t" '{{print $1 "\t" $2 "\t" $2 }}' {input.snp_somatic_hc_filter} | tail -n+2 > {output.snp_somatic_hc_filter_bed}
        """

rule readcount:
    input:
        ref = os.path.join(config["ref_dir"], config["ref_basename"] + ".fa"),
        snp_somatic_hc_filter_bed = "varscan/{subject}.varscan.snp.somatic.hc.filter.bed",
        tumor_bam = lambda wildcards: os.path.join("processed_bams/", config[wildcards.subject]["tumor"] + config["ref_basename"] + ".sorted.chr21.bam")
    output:
        readcounts = "varscan/{subject}.readcounts"
    params:
        bamreadcount = config["bam-readcount"]
    shell:
        """
        {params.bamreadcount} -q 1 -b 20 -f {input.fa} -l {input.snp_somatic_hc_filter_bed} {input.tumor_bam} > {output.readcounts}
        """

rule perl_filter:
    input:
        snp_somatic_hc_filter = "varscan/{subject}.varscan.snp.somatic.hc.filter",
        readcounts = "varscan/{subject}.readcounts"
    output:
        out = "varscan/{subject}.varscan.variants.filter.pass"
    params:
        basename = "varscan/{subject}.varscan.variants.filter",
        perlfilter = config["perl_fp_filter"]
    shell:
        """
        perl {params.perlfilter} {input.snp_somatic_hc_filter} {input.readcounts} --output-basename {params.basename}
        """
