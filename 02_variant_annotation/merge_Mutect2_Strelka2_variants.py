import argparse
from collections import defaultdict

def variant_overlap(gatk, strelka, vep_format):
    """

    :param vcfs:
    :return:
    """
    variants_counter = defaultdict(list)

    # Initialize a summary file to document the number of variants there are in each vcf file
    num_variants_summary = open("number_of_variants_summary.txt", "w")
    header = ["filename", "num_variants"]
    print("\t".join(header), file=num_variants_summary)

    # Initilize an output for vep
    vep_format = open(vep_format, "w")

    # Parse GATK file
    num_variants = 0
    with open(gatk, "r") as f:
        for line in f:
            if not line.startswith("#"):
                num_variants += 1
                cols = line.rstrip("\n").split("\t")
                variant_id = cols[0] + "-" + cols[1] + "-" + cols[3] + "-" + cols[4]
                variants_counter[variant_id].append(1)
        out = ["gatk", str(num_variants)]
        print("\t".join(out), file=num_variants_summary)

    # Parse Strelka file
    num_variants = 0
    with open(strelka, "r") as f: 
        for line in f:
            if not line.startswith("#"):
                num_variants += 1
                cols = line.rstrip("\n").split("\t")
                variant_id = cols[0] + "-" + cols[1] + "-" + cols[3] + "-" + cols[4]
                if len(variants_counter[variant_id]) >= 1:
                    variants_counter[variant_id].append(1)
                else:
                    if len(previous_line) ==0:
                        continue
                    if len(previous_line) > 0:
                        variant_id = cols[0] + "-" + cols[1] + "-" + cols[3] + previous_line[3] + "-" + cols[4] + previous_line[4]
                        if len(variants_counter[variant_id]) >= 1:
                            #print(variant_id)
                            variants_counter[variant_id].append(1)
            previous_line = cols
        out = ["strelka", str(num_variants)]
        print("\t".join(out), file=num_variants_summary)

    at_least_2_vcfs = 0
    for variant_id in variants_counter:
        if len(variants_counter[variant_id]) >= 2:
            cols = variant_id.split("-")
            vep_row = [cols[0], cols[1], ".", cols[2], cols[3], ".", ".", "."]
            print("\t".join(vep_row), file=vep_format)

            at_least_2_vcfs += 1
    out = ["at_least_2_vcfs", str(at_least_2_vcfs)]
    print("\t".join(out), file=num_variants_summary)

def main(args):
    variant_overlap(args.gatk_file, args.strelka_file, args.vep_format_fn)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare input for the program variant effect predictor from the overlap of GATK and Strelka")
    parser.add_argument("--gatk_file", required=True,
                        help="Enter the path to the gatk filename")
    parser.add_argument("--strelka_file", required=True,
                        help= "Enter the path to a combined Strelka indel and snp file")
    parser.add_argument("--vep_format_fn", required=True,
                        help="Enter the path to the output file in VEP format")

    return parser.parse_args()

main(parse_args())
