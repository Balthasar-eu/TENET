#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from jpype import *
import numpy as np
import pathlib
import csv
import argparse
import datetime

# define and parse args

parser = argparse.ArgumentParser(description='Making gene regulatory networks')
parser.add_argument('pairs', help='''The file containing the abtch gene pairs
                    this run of runTE should work on''')
parser.add_argument('output', default="outputs",
                    help='The name of the output file')
parser.add_argument('expression',
                    help='The name of the file containing the expression data')
parser.add_argument('--hist', default=1, type=int,
                    help='History length for transfer entropy calculation')

args = parser.parse_args()

pair_file = args.pairs
output_fname = args.output
historyLength = args.hist
expression_file = args.expression

# load input files
cell_gene_all = np.loadtxt(expression_file)
list_pairs = np.loadtxt(pair_file,dtype=int)

# Change location of jar to match yours:
abspath = pathlib.Path(__file__).parent.absolute()
jarLocation = abspath.joinpath("infodynamics.jar")
# Start the JVM (add the "-Xmx" option with say 1024M if you get crashes due to not enough memory space)
startJVM(getDefaultJVMPath(), "-ea", f"-Djava.class.path={jarLocation}","-Xmx16G")

TEresult=[None] * len(list_pairs)
for num_pair,(p1,p2) in enumerate(list_pairs):
    
    expression_data = [cell_gene_all[p1-1],cell_gene_all[p2-1]]
    
# Create a TE calculator and run it:
    teCalcClass = JPackage("infodynamics.measures.continuous.kernel").TransferEntropyCalculatorKernel
    teCalc = teCalcClass()
    teCalc.setProperty("NORMALISE", "true") # Normalise the individual variables
    teCalc.initialise(historyLength, 0.5) # Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
    resultTemp=[]
    teCalc.setObservations(JArray(JDouble, 1)(expression_data[0]), JArray(JDouble, 1)(expression_data[1]))
    resultTemp.append(teCalc.computeAverageLocalOfObservations())
    teCalc.initialise(historyLength, 0.5) # Use history length 1 (Schreiber k=1), kernel width of 0.5 normalised units
    teCalc.setObservations(JArray(JDouble, 1)(expression_data[1]), JArray(JDouble, 1)(expression_data[0]))
    resultTemp.append(teCalc.computeAverageLocalOfObservations())
    TEresult[num_pair] = [p1,p2] + resultTemp
    if (num_pair % int(len(list_pairs)/2)) == 0:
        print(datetime.datetime.now())

with open(output_fname, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(TEresult)
