#! /usr/bin/python3

import glob

filenames = glob.glob("*_.out")
for i in filenames:
    gene = i[:-5]
    with open(i, "r") as f:
        for line in f:
            if "ls" in line:
                words = line.strip().split()
                length = words[-1]
                print(gene + "\t" + length)


#CODONML (in paml version 4.9h, March 2018)  LOC100643465_cds_aln.fa                                                                                                         Model: several dN/dS ratios for branches,                                                                                                                                   Codon frequency model: F3x4                                                                                                                                                 ns =   2  ls = 408                                                                                                                                                          ^
