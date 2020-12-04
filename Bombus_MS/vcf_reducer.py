#!/usr/bin/python2.7
#remove columns of outliers from vcf

import sys

indices = list(range(0,11)) + [24] + list(range(12,23)) + list(range(26,36))
indices.sort()
#print(indices)
new_words=[]

with open(sys.argv[1],"r") as f:
    for line in f:
        words = line.strip().split()
        new_words = []
        for i in indices:
            new_words.append(words[i])
        new_line = "\t".join(new_words)
        print (new_line)

