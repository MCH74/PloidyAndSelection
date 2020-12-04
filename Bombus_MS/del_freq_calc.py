#! /usr/bin/python3


#use python3 del_freq_calc.py sm sf h uf um

#import sys
import argparse

def del_freq(sm,sf,h,uf,um):
    top = 2 * uf + um
    qd = top / (2 * sf * h)
    qh = top / sm
    qu = top / (2 * sf * h + sm)
    return qd, qh, qu


parser = argparse.ArgumentParser()
parser.add_argument('-sm', action='store', default=0.01,type=float,
                    dest='sm')
parser.add_argument('-sf', action='store',default=0.01,type=float,
                    dest='sf')
parser.add_argument('-dom', action='store',default=0.01,type=float,
                    dest='h')
parser.add_argument('-uf', action='store',default=0.01,type=float,
                    dest='uf')
parser.add_argument('-um', action='store',default=0.01,type=float,
                    dest='um')

args = parser.parse_args()


qd, qh, qu = del_freq(args.sm,args.sf,args.h,args.uf,args.um)
print("diploid frequency: ", round(qd,5))
print("haploid frequency: ", round(qh,5))
print("non-biased frequency: ", round(qu,5))

print("qd:qh ", str(qd/qh))
print("qn:qh ", str(qu/qh))
print("qd:qn ", str(qd/qu))
