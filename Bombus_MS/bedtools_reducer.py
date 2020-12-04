#!/usr/bin/python2.7
#filter bedtools output


with open(sys.argv[1],"r") as f:
    for line in f:
        words = line.strip().split()
        start = words[1]
        end = int(start) + 1
        print (words[0] + "\t" + str(start) + "\t" + str(end) + "\t" + words[2])

