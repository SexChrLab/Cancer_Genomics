import json
import os

data = {}

# fastq path
data["fastq_path"] = "/data/CEM/shared/controlled_access/Beauty/"

# all sample ids
data["all_samples"] = []
with open("/scratch/tphung3/Cancer_Genomics/00_misc/samples_info.csv", "r") as f: #TODO: update the path here
    for line in f:
        if not line.startswith("sampleID"):
            data["all_samples"].append(line.rstrip('\n').split(",")[0])

# populate read group information for each sample
for sample in data["all_samples"]:
    read_group_info = {}
    read_group_info[sample] = {"fq_1": sample + "_1.fastq",
        "fq_2": sample + "_2.fastq",
        "fq_1_fixed": sample + "_1_fixed.fastq",
        "fq_2_fixed": sample + "_2_fixed_fastq",
        "ID": sample,
        "SM": sample,
        "LB": sample,
        "PU": sample,
        "PL": "Illumina"}

    data.update(read_group_info)

# add path to reference
data["ref_dir"] = "/data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome"
data["ref_basename"] = "GRCh38_full_analysis_set_plus_decoy_hla"

with open("somatic_mutation_calling_config.json", "w") as outfile:
    json.dump(data, outfile)

