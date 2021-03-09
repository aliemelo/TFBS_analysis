#!/usr/env/python
import argparse

def main():
    parser = argparse.ArgumentParser(description="Convert FIMO tsv file to gff")
    parser.add_argument("--tsv", action="store", help=".tsv file outputted from FIMO and filtered")
    args = parser.parse_args()

    outfilename = args.tsv + "_converted_to_gff"

    with open(args.tsv, "r") as infile, open(outfilename, "w") as outfile:
        for line in infile:
            if "motif_id" in line:
                continue
            if not (line.startswith("#") or line.startswith("\n")):
                line = line.strip().split("\t")
                pwm = line[0]
                alias = line[1]
                chrom = line[2]
                if "scga7" not in chrom:
                    chrom = "scga7_" + chrom
                start = line[3]
                end = line[4]
                strand = line[5]
                pvalue = line[7]
                qvalue = line[8]
                seq = line[9]

            attribute = "Name={}_{}{};Alias={};ID={}-{}-{};pvalue={};qvalue={};sequence={}".format(pwm, chrom, strand, alias, pwm, alias, chrom, pvalue, qvalue, seq)

            outfile.write(chrom + "\tfimo\tnucleotide_motif\t" + start + "\t" + end + "\t.\t" + strand + "\t.\t" + attribute + "\n" )

if __name__ == "__main__":
    main()
