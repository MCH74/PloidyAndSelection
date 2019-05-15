#! /usr/bin/python3
# -*- coding: utf-8 -*-

import re
import argparse
from collections import defaultdict

def calc_lengths(gff, name):
    summ_dict = defaultdict(list)
    cds_dict = defaultdict(set)
    n_regex = re.escape(name) + "([^,;]+)"

    # reading the gff file and saving cds+exon positions
    with open(gff, 'r') as gf:
        for line in gf:
            if line[0] not in ('\n', '#'):
                llis = line.strip().split()
                classf = llis[2]
                if (classf == 'CDS'):
                    try:
                        gene_ID = re.search(n_regex, llis[8]).group(1)
                        cds_dict[gene_ID].add((int(llis[3]), int(llis[4])))
                    except:
                        print("No matching geneID found for CDS: " + line)

    # calculating gene lengths
    for gene, CDSs in cds_dict.items():
        length = 0
        lCDSs = list(CDSs)
        lCDSs.sort(key=lambda tup: tup[0])
        for cds in lCDSs:
            length += cds[1] - cds[0]
        summ_dict[gene].append(length)

    # calculating intron lengths
    for gene, exons in cds_dict.items():
        introns = []
        intron_len_total = 0
        lexons = list(exons)
        lexons.sort(key=lambda tup: tup[0])
        if (len(lexons) > 1):
            for ind in range(len(lexons)-1):
                introns.append((lexons[ind][1], lexons[ind+1][0]))

            for intron in introns:
                intron_length = (intron[1]-intron[0])-1
                intron_len_total += intron_length

            intron_len_average = float(intron_len_total) / float(len(introns))
            summ_dict[gene].append(intron_len_average)
        else:
            summ_dict[gene].append("N/A")

    return summ_dict

def summary(sum_dict, outfile):
    with open(outfile, 'w') as of:
        of.write("# geneID\tgene length\taverage intron length\n")
        for geneID, rlist in sum_dict.items():
            srlist = [str(v) for v in rlist]
            of.write(geneID + '\t' + '\t'.join(srlist) + '\n')

def main():
    """
    The main function of the program.
    """
    parser = argparse.ArgumentParser(description='Reads a gff file and calculates gene lengths and average intron lengths.')
    parser.add_argument("-g", "--gff", help="The gff file to parse.", type=str, required=True)
    parser.add_argument("-n", "--name", help="The name of the gene identifier field in the gff (default: Parent=).", type=str, default='Parent=')
    parser.add_argument("-o", "--output", help="The name (and path) of the output file in which the results should be written.")
    args = parser.parse_args()

    gene_dict = calc_lengths(args.gff, args.name)

    summary(gene_dict, args.output)


if __name__ == "__main__":
    main()
