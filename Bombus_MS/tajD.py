#!/usr/bin/python2.7
"""
Chrom   pos     Gene    pos     N_ALLELES                                                                                                        NC_015762.1     17946   LOC100649911    603     3       38      -                                                                                NC_015762.1     18132   LOC100649911    417     5       28      - 

CHROM   BIN_START       N_SNPS  TajimaD                                                                                                          NC_015762.1     0       10      -0.374719                                                                                                        NC_015762.1     10000   7       -0.11986                  
"""

#python3 theta.py <freq file> <site counts>
import sys
from collections import defaultdict

taj_dict = {}
with open("tajimasd.Tajima.D","r") as f:
    f.readline()
    for line in f:
        words = line.strip().split()
        chrom = words[0]
        start = int(words[1])
        taj = words[3]
        if(taj == "nan"):
            taj = 0
        else:
            taj = float(taj)
        if (chrom in taj_dict):
            taj_dict[chrom][start] = taj
        else:
            taj_dict[chrom] = {}
            taj_dict[chrom][start] = taj

mean_taj = defaultdict  (list)
with open(sys.argv[1],"r") as f:
    f.readline()
    for line in f:
        words = line.strip().split()
        if (int(words[4]) > 0 and int(words[4]) != int(words[5]) ):
            chrom = words[0]
            base = int(words[1])
            gene = words[2]
            start = int(base / 10000) * 10000
            #print("chrom:", chrom, "start:", start)
            mean_taj[gene].append(taj_dict[chrom][start])
            #if(gene in mean_taj):
             #   mean_taj[gene].append(taj_dict[chrom][start])
            #else:
             #   mean_taj[gene] = [taj_dict[chrom][start]]
              #  print(mean_taj)
lengths = {}
with open(sys.argv[2],"r") as f:
    for line in f:
        words = line.strip().split()
        lengths[words[0]] = words[1]

for gene in mean_taj:
    if(gene in lengths):
        taj = sum(mean_taj[gene]) / int(lengths[gene])
        print(gene + "\t" + str(taj) )

