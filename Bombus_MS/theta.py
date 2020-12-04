#!/usr/bin/python2.7
"""
Chrom pos Gene pos      N_ALLELES   Strand
NC_015762.1    17946    LOC100649911 603 3  38  -
NC_015762.1    18132    LOC100649911 417 5  28  -
"""

#python3 theta.py <freq file> <site counts>
import sys

harm = {}
with open("harmonic_series","r") as f:
    for line in f:
        words = line.strip().split()
        harm[words[0]] = words[1]

theta_dict = {}
with open(sys.argv[1],"r") as f:
    f.readline()
    for line in f:
        words = line.strip().split()
        gene = words[2]
        alleles = int(words[4])
        total = int(words[5])
        if(alleles < total and alleles > 0):
            theta = 1 / float(harm[str(total - 1)])
        else:
            theta = 0
        if(gene in theta_dict):
            theta_dict[gene] = theta + theta_dict[gene]
        else:
            theta_dict[gene] = theta
        
lengths = {}
with open(sys.argv[2],"r") as f:
    for line in f:
        words = line.strip().split()
        lengths[words[0]] = words[1]

for gene in theta_dict:
    if(gene in lengths):
        theta = theta_dict[gene] / int(lengths[gene])
        print(gene + "\t" + str(theta) )

