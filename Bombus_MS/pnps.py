#!/usr/bin/python2.7
#calculate pN/pS from allele counts file
"""
Chrom pos Gene pos      N_ALLELES    Strand
NC_015762.1    17946    LOC100649911 603 3  38  -
NC_015762.1    62989    LOC100648154 729 2  42  +
"""

#python3 pnps.py 
import sys

def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 

###number of sites
four_counts = {}
zero_counts = {}
gene_lengths = {}

with open("four_fold_site_counts","r") as f:
    for line in f:
        words = line.strip().split()
        four_counts[words[0]] = words[1]

with open("zero_fold_site_counts","r") as f:
    for line in f:
        words = line.strip().split()
        zero_counts[words[0]] = words[1]

with open("gene_lengths","r") as f:
    for line in f:
        words = line.strip().split()
        gene_lengths[words[0]] = words[1]

##get genes shared by all
shared_genes = intersection(  list( four_counts.keys() ) , list( gene_lengths.keys() ) )
#shared_genes = intersection( shared_genes, list( gene_lengths.keys() ) )
#print(shared_genes)

###number of segregating sites
alleles4 = {}
alleles0 = {}
alleles_tot = {}

with open("allele_counts_cds.4","r") as f:
    f.readline()
    for line in f:
        words = line.strip().split()
        if(int(words[4]) > 0 and int(words[4]) < int(words[5])):
            if( words[2] in alleles4):
                alleles4[words[2]] = alleles4[words[2]] + 1
            else:
                alleles4[words[2]] = 1

with open("allele_counts_cds.0","r") as f:
    f.readline()
    for line in f:
        words = line.strip().split()
        if(int(words[4]) > 0 and int(words[4]) < int(words[5])):
            if( words[2] in alleles0):
                alleles0[words[2]] = alleles0[words[2]] + 1
            else:
                alleles0[words[2]] = 1

with open("allele_counts_cds.nonredundant","r") as f:
    f.readline()
    for line in f:
        words = line.strip().split()
        if(int(words[4]) > 0 and int(words[4]) < int(words[5])):
            if( words[2] in alleles_tot):
                alleles_tot[words[2]] = alleles_tot[words[2]] + 1
            else:
                alleles_tot[words[2]] = 1


print("Gene\tpT\tpN\tpS\tpNpS")

for gene in shared_genes:
    if(gene in alleles_tot):
        pt = int(alleles_tot[gene]) / int(gene_lengths[gene])
    else:
        pt = 0
    if(gene in alleles0):
        pn = int(alleles0[gene]) / int(zero_counts[gene])
    else:
        pn = 0
    if(gene in alleles4):
        ps = int(alleles4[gene]) / int(four_counts[gene])
    else:
        ps = 0
    if(ps == 0):
        pnps = "NA"
    else:
        pnps = str(pn/ps)
    print(gene + "\t" + str(pt) + "\t" + str(pn) + "\t" + str(ps) + "\t" + pnps)
    #print(gene + "\t" + pnps)

