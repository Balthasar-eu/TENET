#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import os
import pathlib
import tempfile

###############################################################################
############################ parameter  collection ############################
###############################################################################

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Making gene regulatory networks')
    parser.add_argument('inputcsv', help='''The csv file containing the reads with
                        genes in rows and samples columns''')
    parser.add_argument('trajectory_file',
                        help='The trajectory_file containing RNA velocity information')
    parser.add_argument('cell_select',
                        help='The file containing the selected cells')
    parser.add_argument('output', default="outputs",
                        help='The directory in which the output should be written')
    parser.add_argument('--transposed', default=0, type=bool,
                        help='''When the input csv is transposed and has genes in
                        rows and samples in columns set this to one''')
    parser.add_argument('--jobs', default=1, type=int,
                        help='How many jobs should be run in parallel')
    parser.add_argument('--hist', default=1, type=int,
                        help='History length for mpirun')
    #parser.add_argument('--tf', default=0, type=bool,
    #                    help='Run with transcription factors only')
    
    args = parser.parse_args()
    
    inputcsv = args.inputcsv
    trajectory_file = args.trajectory_file
    cell_select = args.cell_select
    output_dir = args.output
    JOBNUM = args.jobs
    HIST_LENGTH = args.hist
    #TF = args.tf
    TRANSPOSED = args.transposed
    
###############################################################################
################################ preprocessing ################################
###############################################################################

branch=np.loadtxt(cell_select,dtype=int)
usecols = [i+1 for i,b in enumerate(branch == 1) if b]

if TRANSPOSED:
    with open(inputcsv) as f:
        #skip first empty field
        gene_names = f.readline().strip().split(",")[1:]
    
    expression_data=np.loadtxt(inputcsv,skiprows=1,delimiter=",")[usecols]
    
else:
    gene_names = []
    with open(inputcsv) as f:
        #skip header
        f.readline()
        for l in f:
            row = l.strip().split(",")
            gene_names.append(row[0])
    
    expression_data=np.loadtxt(inputcsv,
                               dtype=float,
                               usecols=usecols,
                               skiprows=1,
                               delimiter=",").T

num_gene=expression_data.shape[0]

r = np.arange(1,num_gene+1)
x1,x2 = np.meshgrid(r,r)
i1,i2 = np.triu_indices(num_gene,k=1)
indx = np.array([x2[i1,i2],x1[i1,i2]]).T

pairs = np.array_split(indx,JOBNUM)

###############################################################################
#################################### TENET ####################################
###############################################################################

trajectory=np.loadtxt(trajectory_file)
trajectory=trajectory[branch==1]
trajectorySortIndex=np.argsort(trajectory)

expression_data = expression_data[trajectorySortIndex]
#np.savetxt("cell_gene_trsps.tsv",expression_data)

pairnames = []
for i, pair in enumerate(pairs):
    with tempfile.NamedTemporaryFile(mode="w+b", delete=False) as temp:
        pairnames.append(temp.name)
        np.savetxt(temp,pair,fmt="%d")

pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

cmd = 'time mpirun' 

executable = pathlib.Path(__file__).parent.absolute().joinpath("runTE_wrapped.py")

with tempfile.NamedTemporaryFile(mode="w+b") as ntf:
    np.savetxt(ntf,expression_data)
    ntf.flush()
    
    outnames = []
    for i,pn in enumerate(pairnames):
        out = str(pathlib.Path(output_dir).absolute().joinpath(f"TE_out_{i:02d}.csv"))
        cmd += f" -np 1 {executable} {pn} {out} {ntf.name} --hist {HIST_LENGTH} :"
        outnames.append(out)

    print(cmd,"\n")

    # RUN THE COMMANDS
    os.system(cmd)

TEmatrix = np.zeros((num_gene,num_gene))

for fname in outnames:
    ifile = open(fname)
    for line in ifile:
        temp=line.strip().split(",")
        if len(temp)>3:
            TEmatrix[int(temp[1])-1][int(temp[0])-1]=float(temp[3])
        else:
            TEmatrix[int(temp[0])-1][int(temp[1])-1]=float(temp[2])

out = str(pathlib.Path(output_dir).joinpath("TE_result_matrix.txt"))
with open(out,"w") as ofile:
    ofile.write("TE")
    for g in gene_names:
        ofile.write("\t"+g)
    for i,g in enumerate(gene_names):
        ofile.write("\n"+g)
        for j in range(num_gene):
            ofile.write("\t"+str(TEmatrix[i][j]))    

