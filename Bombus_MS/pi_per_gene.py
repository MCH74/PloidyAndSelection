#!/usr/bin/python2.7
"""
Chrom   pos     Gene    pos     pi                                                                                                               NC_015762.1     17946   LOC100649911    603     0.14936                                                                                          NC_015762.1     18132   LOC100649911    417     0.304233                                                                                         N"""

#python3 <>.py pi file> <site counts>
import sys

pi_dict = {}
with open(sys.argv[1],"r") as f:
    f.readline()
    for line in f:
        words = line.strip().split()
        gene = words[2]
        if (words[4] == "-nan"):
            pi = 0
        else:
            pi = float(words[4])
        if(gene in pi_dict):
            pi_dict[gene] = pi + pi_dict[gene]
        else:
            pi_dict[gene] = pi
        
lengths = {}
with open(sys.argv[2],"r") as f:
    for line in f:
        words = line.strip().split()
        lengths[words[0]] = words[1]

for gene in pi_dict:
    if(gene in lengths):
        pi = pi_dict[gene] / int(lengths[gene])
        print(gene + "\t" + str(pi) )

