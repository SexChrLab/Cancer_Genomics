import json
from collections import defaultdict
import os

data = {}

# fastq path
data["fastq_path"] = "/data/CEM/shared/controlled_access/Beauty/"

# all sample ids
data["all_samples"] = []

all_subjects = set()
data["all_subjects"] = []

subjectID_normal_tumor_info = defaultdict(list)
with open("/scratch/tphung3/Cancer_Genomics/00_misc/samples_info.csv", "r") as f: #TODO: update the path here
    for line in f:
        if not line.startswith("subjectID"):
            items = line.rstrip("\n").split(",")
            data["all_samples"].append(items[1])
            all_subjects.add(items[0])

            subjectID_normal_tumor_info[items[0]].append(items[1]) #Assuming that normal is before tumor

for i in all_subjects:
    data["all_subjects"].append(i)

# add information about normal and tumor
for i in subjectID_normal_tumor_info:
    info = {}
    info[i] = {"normal": subjectID_normal_tumor_info[i][0],
         "tumor": subjectID_normal_tumor_info[i][1]}
    data.update(info)

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

# add path to reference
data["ref_dir"] = "/data/CEM/shared/public_data/references/1000genomes_GRCh38_reference_genome"
data["ref_basename"] = "GRCh38_full_analysis_set_plus_decoy_hla"

# add path to software
data["varscan_path"] = "/home/tphung3/softwares/VarScan.v2.3.9.jar"
data["gatk_path"] = "/home/tphung3/softwares/gatk-4.1.8.1/gatk"
data["bam-readcount"] = "/home/tphung3/softwares/bam-readcount"
data["perl_fp_filter"] = "/home/tphung3/softwares/fpfilter-2.pl"

with open("somatic_mutation_calling_config.json", "w") as outfile:
    json.dump(data, outfile)

