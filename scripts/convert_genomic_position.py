#!/usr/bin/python
import argparse
import re
from collections import defaultdict


def create_ref_dict(bedfile):
    out = defaultdict(list)
    with open(bedfile, 'r') as infile:
        for lines in infile:
            line = lines.strip().split()
            gene = re.match(r"ID=(.*);Name", line[3])
            gene = gene.group(1)
            if "mikado" in gene:
                gene = gene[7:]
            if "evm" in gene:
                gene = gene[7:]
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            strand = line[5]
            out[gene].append(chrom)
            out[gene].append(start)
            out[gene].append(end)
            out[gene].append(strand)
    return out


def main():
    parser = argparse.ArgumentParser(description="Script to change the feature coordinates from sequence-based to chromosome-based")
    parser.add_argument("-f", action="store", help="GFF file to be converted", required=True)
    parser.add_argument("--bed", action="store", help="BED file with the promoter region coordinates")
    args = parser.parse_args()

    # create dictionary with reference positions from bed file
    ref_dict = create_ref_dict(args.bed)

    # convert positions and write to outfile
    with open(args.f, 'r') as infile, open(args.f + '.converted_coords', 'w') as outfile:
        for line in infile:
            line = line.strip().split()
            if "ID=" in line[0]:
                gene = line[0][3:]
            else:
                gene = line[0]
            start = int(line[3])
            end = int(line[4])
            strand = line[6]
            info = line[8]

            if gene not in ref_dict:
                continue

            chrom, ref_start, ref_end, ref_strand = [ref_dict[gene][i] for i in range(4)]

            new_start = str(ref_start + start)
            new_end = str(ref_start + end)

            outfile.write("{}\tfimo\tnucleotide_motif\t{}\t{}\t.\t{}\t.\t{}\n".format(chrom,
            new_start, new_end, strand, info))

if __name__ == "__main__":
    main()
