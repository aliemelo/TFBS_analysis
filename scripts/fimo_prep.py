#!/usr/env/python
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--file", action="store", help="")
    parser.add_argument("--alt", action="store")
    args = parser.parse_args()

    tfbs_dict = defaultdict(list)

    with open(args.file, 'r') as infile:
        for line in infile:
            line = line.strip().split(",")
            seq = line[0]
            motif = line[5]
            if "/" in motif:
                motif = motif.split("/")
                motif = set(motif)
                motif = '/'.join(sorted(motif))
            tfbs_dict[seq].append(motif)

    # TSV com sequencia e sítios COM alternativa
    if args.alt:
        with open('resultados_mapeamento.tsv', 'w') as ofile:
            for seq, motifs in tfbs_dict.items():
                motifs = ','.join(motifs)
                ofile.write(seq + "\t" + motifs + "\n")
    else:
        # TSV com sequencia e sítios SEM alternativa
        with open('resultados_mapeamento.tsv', 'w') as out_file:
            for seq, motifs in tfbs_dict.items():
                motifs = '/'.join(motifs)
                motifs = motifs.split("/")
                motifs = list(set(motifs))
                motifs = ','.join(motifs)
                out_file.write(seq + "\t" + motifs + "\n")


if __name__ == "__main__":
    main()