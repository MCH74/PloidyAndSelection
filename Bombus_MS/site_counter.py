#!/usr/bin/python2.7
#count number of 4- and 0-fold sites per gene
import sys


with open(sys.argv[1],"r") as f:
    for line in f:
        words = line.strip().split()
        gene = words[0]
        sites = len(words[1:])
        print (gene + "\t" + str(sites))

