#!/usr/env/python
import argparse


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", "--file", action="store", help="", required=True)
    args = parser.parse_args()

    output = args.file + ".start.bed"
    with open(args.file, "r") as gff_file, open(output, 'w') as outfile:
        for line in gff_file:
            columns = line.strip().split("\t")
            seqname = columns[0]
            start = columns[3]
            end = columns[4]
            score = columns[5]
            strand = columns[6]
            attrb = columns[8]

            if strand == '+':
                end = str(int(start) + 2)
                outfile.write(seqname + "\t" + start + "\t" + end + "\t" + attrb + "\t" + score + "\t" + strand + "\n")
            else:
                start = str(int(end) - 2)
                outfile.write(seqname + "\t" + start + "\t" + end + "\t" + attrb + "\t" + score + "\t" + strand + "\n")


if __name__ == "__main__":
    main()