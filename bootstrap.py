#!/usr/bin/python3
import os
import numpy as np, scipy.stats as st

def conf(a): #calculate confidence interval with array a as argument
    conf_int = st.t.interval(0.95, len(a)-1, loc=np.mean(a), scale=st.sem(a))
    return conf_int

files = os.listdir("Bootstrapping/")
sos1 = []
sos2 = []
sos3 = []
sos4 = []
sos5 = []
gammaz = []
fww = []
basic = []
for element in files:
    infile = "Bootstrapping/" + element
#    print(infile)
    with open(infile, "r") as f:
        words = f.read().split()
    n = -1
    for word in words:
        n += 1
        if "likelihood" in word:
            sos1.append(float(words[int(n)+15][12:]))
            sos2.append(float(words[int(n)+16][11:]))
            sos3.append(float(words[int(n)+17][10:]))
            sos4.append(float(words[int(n)+18][9:]))
            sos5.append(float(words[int(n)+19][7:]))
            basic.append(float(words[int(n)+27][8:-2]))
            fww.append(float(words[int(n)+28][12:-2]))
            gammaz.append(float(words[int(n)+29][12:]))

print(np.mean(sos1), conf(sos1)[0], conf(sos1)[1],)
print(np.mean(sos2), conf(sos2)[0], conf(sos2)[1]) 
print(np.mean(sos3), conf(sos3)[0], conf(sos3)[1]) 
print(np.mean(sos4), conf(sos4)[0], conf(sos4)[1]) 
print(np.mean(sos5), conf(sos5)[0], conf(sos5)[1]) 
print(np.mean(basic), conf(basic)[0], conf(basic)[1]) 
print(np.mean(fww), conf(fww)[0], conf(fww)[1]) 
print(np.mean(gammaz), conf(gammaz)[0], conf(gammaz)[1]) 

"""         print(n)
            print(words[int(n)+15][12:])
            print(words[int(n)+16][11:])
            print(words[int(n)+17][10:])
            print(words[int(n)+18][9:])
            print(words[int(n)+19][7:])
            print(words[int(n)+27][8:-2])
            print(words[int(n)+28][12:-2])
            print(words[int(n)+29][12:])"""
