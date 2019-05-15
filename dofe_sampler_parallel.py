#!/usr/bin/python3
import numpy as np
import pandas as pd
import subprocess
import sys

"""
take full dofe input table and generate one normal output from sum of all genes and n random samples with replacement for bootstrapping
use:
python3 dofe_sample.py infile groupname number_of_sample suffix
"""

#generate random sample of indices for bootstrapping
def sample(length): #argument is length of table, i.e. number of rows
    indices = np.random.choice(list(range(0,length)), size = length, replace=True)
    return indices

#run dofe command with grapes on input file
def dofe_cmd(dofe_infile,group):
    proc = subprocess.Popen(['/global/group/home/harrisom/Leicester/grapes-master/grapes/grapes', '-in', dofe_infile, '-out', 'out.txt', '-model', 'GammaZero'], stderr = subprocess.PIPE, stdout = subprocess.PIPE)
    stdout,stderr = proc.communicate()
    return proc.returncode, stdout, stderr

infile = sys.argv[1]
group = sys.argv[2]
n = int(sys.argv[3])
suffix = sys.argv[4]

dofe_df = pd.read_csv(infile,sep="\t",header=None)

#first standard output
dofe_in = dofe_df.iloc[:,2:].sum().values.tolist()
s = [str(i ) for i in dofe_in]
dofe_in = '\t'.join(s)
dofe_infile = "dofe_in_total.txt"
with open(dofe_infile,"w") as f:
    f.write(group + "\n" + "summed" + "\t" + "80" + "\t" + str(dofe_in)) 
code, out, err = dofe_cmd(dofe_infile,group)
outfile = group + "_dofe.out"
with open(outfile,"w") as f:
    f.write("out: '{}'".format(out))

"""
length = dofe_df.shape[0]
for i in range(1,n+1):
    outfile = "Bootstrapping/" + group + "_BS_" + str(i) + "_" + suffix + ".txt"
    indices = sample(length)
    dofe_in = dofe_df.iloc[indices,2:].sum().values.tolist()
    s = [str(i ) for i in dofe_in]
    dofe_in = '\t'.join(s)
    dofe_infile = "dofe_in_" + str(i) + "_" + suffix + ".txt"
    with open(dofe_infile,"w") as f:
        f.write(group + "\n" + "summed" + "\t" + "80" + "\t" + str(dofe_in)) 
    code, out, err = dofe_cmd(dofe_infile,group)
    with open(outfile,"w") as f:
        f.write("out: '{}'".format(out))



"""
