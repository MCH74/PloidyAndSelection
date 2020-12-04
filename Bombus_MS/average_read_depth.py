#! /usr/bin/python3
#SCAFF   POS     READS                                                                                                                                                       NC_015762.1     389     74                                                                                                                                                  NC_015762.1     394     74                                                                                                                                                  NC_015762.1     544     104                                                                                                                                                 NC_015762.1     664     202                                                                                                                                                 
#harrisom@ebbrome:~/Leicester/Bt_Genome/Star/Poly2/2pass/New/WO_Outliers$ head male_gene_coordinates                                                                         NC_015762.1     279470  285561                                                                                                                                              NC_015762.1     365916  368227                                                                                                                                              NC_015762.1     385861  388241                                                                                                                                             

#use: python3 script <read lengths> <gene list>

import sys

coordinates = {}
with open(sys.argv[2], "r") as f:
    for line in f:
        words = line.strip().split()
        scaff = words[0]
        start = int(words[1])
        end = int(words[2]) + 1
        if scaff not in coordinates:
            coordinates[scaff] = []
        for i in range(start, end):
            coordinates[scaff].append(i)
runsum = 0
i = 0
with open(sys.argv[1], "r") as f:
    f.readline()
    for line in f:
        words = line.strip().split()
        scaff = words[0]
        pos = int(words[1])
        length = int(words[2])
        if scaff in coordinates and pos in coordinates[scaff]:
            i = i + 1
            runsum = runsum + length
print(str(runsum/i))
