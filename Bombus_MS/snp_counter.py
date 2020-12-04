#!/usr/bin/python2.7
#count snps per gene

#python3 theta.py <freq file> <genelist>
import sys

snp_dict = {}
with open(sys.argv[2],"r") as f:
    for line in f:
        gene = line.strip()
        snp_dict[gene] = 0
with open(sys.argv[1],"r") as f:
   # f.readline()
    for line in f:
        words = line.strip().split()
        gene = words[2]
        alleles = int(words[4])
        total = int(words[5])
        if(alleles < total and alleles > 0):
            if(gene in snp_dict):
                snp_dict[gene] = snp_dict[gene] + 1
            else:
                snp_dict[gene] = 1
        
for gene in snp_dict:
    print(gene + "\t" + str(snp_dict[gene]) )

