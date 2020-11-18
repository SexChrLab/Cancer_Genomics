import argparse
from collections import defaultdict

def select_variants_at_least_two_callers(vcfs, vep_format):
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

    for vcf in vcfs:
        num_variants = 0
        with open(vcf, "r") as f:
            for line in f:
                if not line.startswith("#"):
                    num_variants += 1
                    cols = line.rstrip("\n").split("\t")
                    variant_id = cols[0] + "-" + cols[1] + "-" + cols[3] + "-" + cols[4]
                    variants_counter[variant_id].append(1)
        out = [vcf, str(num_variants)]
        print("\t".join(out), file=num_variants_summary)
    at_least_2_vcfs = 0
    for variant_id in variants_counter:
        if len(variants_counter[variant_id]) >= 2:
            cols = variant_id.split("-")
            vep_row = [cols[0][3:], cols[1], ".", cols[2], cols[3], ".", ".", "."]
            print("\t".join(vep_row), file=vep_format)

            at_least_2_vcfs += 1
    out = ["at_least_2_vcfs", str(at_least_2_vcfs)]
    print("\t".join(out), file=num_variants_summary)

def convert_to_vep_format_one_caller(vcf, vep_format):
    vep_format = open(vep_format, "w")
    with open(vcf, "r") as f:
        for line in f:
            if not line.startswith("#"):
                cols = line.rstrip("\n").split("\t")
                out = [cols[0][3:], cols[1], ".", cols[3], cols[4], ".", ".", "."]
                print("\t".join(out), file=vep_format)

def main(args):
    vcfs_list = args.vcf_filenames.split(',')
    if len(vcfs_list) > 1:
        select_variants_at_least_two_callers(vcfs_list, args.vep_format_fn)
    else:
        convert_to_vep_format_one_caller(vcfs_list[0], args.vep_format_fn)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare input for the program variant effect predictor from one or more vcf files")
    parser.add_argument("--vcf_filenames", required=True,
                        help="Enter the path to one vcf file or paths to many vcf files separated by commas")
    parser.add_argument("--vep_format_fn", required=True,
                        help="Enter the path to the output file in VEP format")

    return parser.parse_args()

main(parse_args())