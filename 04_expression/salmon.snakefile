#! importing join
from os.path import join

# Workflow for quasi-quantification of RNAseq read counts with Salmon in non-alignment-based mode.

# Configuration file
configfile: "somatic_mutation_calling_config.json"

# Tools
SALMON = "salmon"

# Reference genome files: XX with Y chromosome masked, XY with both Y-chromosomal PAR masked
SALMON_INDEX = "/data/CEM/shared/public_data/references/ensemble_GRCh37.101/Homo_sapiens.GRCh37.cdna.all"

# Directories
FQ_DIR = "/scratch/eknodel/Cancer_Genomics/01_somatic_mutation_calling/trimmed_fastqs/" # path to directory with trimmed FASTQ files
SALMON_DIR = "/scratch/eknodel/Cancer_Genomics/04_expression/" # path to directory for Salmon output files

# Samples
SAMPLES = config["RNA"]

rule all:
    input:
        # Defining the files that snakemake will attempt to produce as an output.
        # If there is no rule defined to produce the file, or if the file already
        # exists, snakemake will throw "Nothing to be done"
        expand(SALMON_DIR + "{sample}_salmon_quant/", SALMON_DIR=SALMON_DIR, sample=SAMPLES),

rule salmon_quant_paired:
    input:
        R1 = os.path.join(FQ_DIR, "{sample}_trimmomatic_trimmed_paired_1.fastq"),
        R2 = os.path.join(FQ_DIR, "{sample}_trimmomatic_trimmed_paired_2.fastq")
    output:
        OUTPUT = os.path.join(SALMON_DIR, "{sample}_salmon_quant/"),
    params:
        SALMON = SALMON,
        SALMON_INDEX = SALMON_INDEX,
        LIBTYPE = "A", # LIBTYPE A for automatic detection of library type
        threads = 8
    message: "Quantifying {wildcards.sample} transcripts with Salmon."
    run:
        shell("{params.SALMON} quant -i {params.SALMON_INDEX} -l {params.LIBTYPE} -1 {input.R1} -2 {input.R2} -o {output.OUTPUT}")

