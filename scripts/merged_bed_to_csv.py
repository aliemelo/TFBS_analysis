#!/usr/env/python
import argparse
import re


def main():
    parser = argparse.ArgumentParser(description="Create .csv file with information about the mapping of JASPAR database")
    parser.add_argument("--bed", action="store", help="Merged bed")
    args = parser.parse_args()

    with open(args.bed, 'r') as infile, open(args.bed+'.csv', 'w') as outfile, open("fimo_results_excluded_promoters.txt", "w") as excl_file:
        for line in infile:
            line = line.strip().split("\t")
            promoter = line[0]
            if "," in promoter:
                excl_file.write(promoter + "\n")
                continue
            start = line[1]
            end = line[2]
            strand = line[3]
            attbs = line[4].split(",")
            motifs = []
            seqs = []
            aliases = []
            for item in attbs:
                motif = re.match(r'Name=(.*)_scga7.*;Alias', item).group(1)
                motifs.append(motif)
                seq = re.match(r'.*sequence=(.*)', item).group(1)
                seqs.append(seq)
                alias = re.match(r'.*Alias=(.*);ID', item).group(1)
                aliases.append(alias)
            motifs.sort()
            aliases.sort()
            m = '/'.join(motifs)
            s = '/'.join(seqs)
            a = '/'.join(aliases)
            outfile.write(promoter + "," + start + "," + end + "," + strand + "," + m + "," + a + "," + s + "\n")


if __name__ == "__main__":
    main()