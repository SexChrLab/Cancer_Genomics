# Setting up filenames here:
from os.path import join
#configfile: "full_hla_types.json"
configfile: "somatic_mutation_calling_config.json"

rule all:
    input:
        expand("netMHCpan/{sample}_9_netmhc.xsl", sample=config["all_subjects"]), # Calculate dissociation constant
        expand("netMHCstabpan/{sample}_9_netmhcstab.xsl", sample=config["all_subjects"]) # Calculate binding stability

rule prepare_input: 
    input: 
        peptides = os.path.join("../02_variant_annotation/peptides/{sample}_vep.17.peptide")
    output: 
        peptides = os.path.join("../02_variant_annotation/peptides/{sample}.peptide_formatted")
    shell:
        """
        cat {input.peptides} | sed '/>WT/,+1 d' > {output.peptides}
        """

rule netMHC:
    input:
        peptides = os.path.join("../02_variant_annotation/peptides/{sample}.peptide_formatted")
    output:
        netMHC_9  = os.path.join("netMHCpan/{sample}_9_netmhc.xsl"),
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"][0]
    shell:
        """
        netMHCpan -a {params.hla} -f {input.peptides} -BA -s -xls -l 9  -xlsfile {output.netMHC_9};
        sed -i  '1d' {output.netMHC_9}
        """

rule netMHCstab:
    input:
        peptides = os.path.join("../02_variant_annotation/peptides/{sample}.peptide_formatted")
    output:
        netMHCstab_9  = os.path.join("netMHCstabpan/{sample}_9_netmhcstab_polysolver.xsl"),
    params:
        hla = lambda wildcards: config[wildcards.sample]["hla"][0]
    shell:
        """
        netMHCstabpan -a {params.hla} -f {input.peptides} -s -xls -l 9  -xlsfile {output.netMHCstab_9};
        """
