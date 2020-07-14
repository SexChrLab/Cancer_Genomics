import json
import os

data = {}

# fastq path
data["fastq_path"] = "/data/CEM/shared/controlled_access/Beauty/"

# all sample ids
data["all_samples"] = []
with open("samples_info.csv", "r") as f:
    for line in f:
        if not line.startswith("sampleID"):
            data["all_samples"].append(line.rstrip('\n').split(",")[0])

# populate read group information for each sample
for sample in data["all_samples"]:
    read_group_info = {}
    read_group_info[sample] = {"fq_1": sample + "_1.fastq",
        "fq_2": sample + "_2.fastq",
        "ID": sample,
        "SM": sample,
        "LB": sample,
        "PU": sample,
        "PL": "Illumina"}

    data.update(read_group_info)

with open("somatic_mutation_calling_config.json", "w") as outfile:
    json.dump(data, outfile)

