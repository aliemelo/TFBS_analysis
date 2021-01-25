#!/usr/bin/env python3

import argparse


def main():
    parser = argparse.ArgumentParser(description="Convert background files from MotifSuite format to MEME format")
    parser.add_argument("--file", action="store", help="Background model in MotifSuite format")
    args = parser.parse_args()

    with open(args.file, "r") as infile, open(args.file+"_meme_format", "w") as outfile:
        outfile.write("#\torder 0\n")
        for line in infile:
            if line.startswith("#snf"):
                a, c, g, t = next(infile).split()
        outfile.write("A\t{}\n".format(a))
        outfile.write("C\t{}\n".format(c))
        outfile.write("G\t{}\n".format(g))
        outfile.write("T\t{}\n".format(t))


if __name__ == "__main__":
    main()
            
