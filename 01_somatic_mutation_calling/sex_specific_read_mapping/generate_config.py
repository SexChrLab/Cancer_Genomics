import json
from collections import defaultdict
import os
import argparse

parser = argparse.ArgumentParser(description="Generate config files")
parser.add_argument("--fastq_path",required=True,help="Path to where all the fastqs are. For example: /data/CEM/shared/controlled_access/Beauty/")
parser.add_argument("--sample_info",required=True,help="Path to the sample info. For example: /scratch/tphung3/Cancer_Genomics/00_misc/samples_info.csv")
parser.add_argument("--female_ref_dir",required=True,help="Path to the directory where the female references are. For example: /data/CEM/shared/public_data/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC")
parser.add_argument("--female_ref_basename",required=True,help="Input the basename of the female reference. For example: GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded")
parser.add_argument("--male_ref_dir",required=True,help="Path to the directory where the male references are. For example: /data/CEM/shared/public_data/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC")
parser.add_argument("--male_ref_basename",required=True,help="Input the basename of the male reference. For example: GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded")
parser.add_argument("--ref_basename",required=True,help="Input a basename that is sex neutral. For example: GCA_009914755.4_CHM13_T2T_v2.0")
parser.add_argument("--varscan_path",required=True,help="Path to varscan. For example: /home/tphung3/softwares/VarScan.v2.3.9.jar")
parser.add_argument("--gatk_path",required=True,help="Path to GATK. For example: /home/tphung3/softwares/gatk-4.1.8.1/gatk")
parser.add_argument("--strelka",required=True,help="Path to strelka. For example: /home/tphung3/softwares/miniconda3/envs/cancer/share/strelka-2.9.10-0/bin/configureStrelkaSomaticWorkflow.py")
parser.add_argument("--bam_readcount",required=True,help="Path to bam-readcount. For example: /home/tphung3/softwares/bam-readcount")
parser.add_argument("--perl_fp_filter",required=True,help="Path to perl fp filter. For example: /home/tphung3/softwares/fpfilter-2.pl")

args = parser.parse_args()

data = {}

# fastq path
data["fastq_path"] = args.fastq_path

# all sample ids
data["all_samples"] = []
data["DNA"] = []
data["RNA"] = []
data["female"] = []
data["male"] = []

all_subjects = set()
data["all_subjects"] = []

subjectID_normal_tumor_info = defaultdict(list)
with open(args.sample_info, "r") as f: #TODO: update the path here
    for line in f:
        if not line.startswith("subjectID"):
            items = line.rstrip("\n").split(",")
            if items[4] == "DNA":
                data["DNA"].append(items[1])
                all_subjects.add(items[0])
                if items[5] == "male":
                    data["male"].append(items[1])
                if items[5] == "female":
                    data["female"].append(items[1])

                subjectID_normal_tumor_info[items[0]].append(items[1]) #Assuming that normal is before tumor
            elif items[4] == "RNA":
                data["RNA"].append(items[1])
            data["all_samples"].append(items[1])

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
    read_group_info[sample] = {"fq_1": sample + "_1.fastq.gz",
        "fq_2": sample + "_2.fastq.gz",
        "ID": sample,
        "SM": sample,
        "LB": sample,
        "PU": sample,
        "PL": "Illumina"}

    data.update(read_group_info)

# add path to reference
data["female_ref_dir"] = args.female_ref_dir
data["female_ref_basename"] = args.female_ref_basename
data["male_ref_dir"] = args.male_ref_dir
data["male_ref_basename"] = args.male_ref_basename
data["ref_basename"] = args.ref_basename

# add path to software
data["varscan_path"] = args.varscan_path
data["gatk_path"] = args.gatk_path
data["strelka"] = args.strelka
data["bam-readcount"] = args.bam_readcount
data["perl_fp_filter"] = args.perl_fp_filter

with open("somatic_mutation_calling_config.json", "w") as outfile:
    json.dump(data, outfile)

