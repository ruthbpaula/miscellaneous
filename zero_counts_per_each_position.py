# -*- coding: utf-8 -*-
"""
Created on Thu May 27 10:45:12 2021

@author: Ruth De Paula
"""

import pandas as pd
import numpy as np
import math
import numpy.matlib
import sys

# # arg1 = sys.argv[1]
arg1 = '500_second_splice_g4-22_toPlot_2.txt' # input file with counts

file_toPlot = pd.read_csv(arg1, delimiter = '\t', names = ['coord', 'counts']) # ficticious header to access columns with (.column_name)
file_toPlot_numpy = file_toPlot.to_numpy()

# we need a first comparison to start the loop
index = 0
output_dataframe = np.empty([0,1])

sys.stdout = open("500_second_splice_g4-22_toPlot_2_output.txt", "w")

for i in range(-5000,5001):
    
    if (i == file_toPlot_numpy[index,0]):
#        row = np.concatenate(i, "\t", file_toPlot_numpy[index,1])
#        output_dataframe = np.vstack([output_dataframe, row])
        print(i, "\t", file_toPlot_numpy[index,1])
        index += 1
    
    else:
#        row = np.concatenate(i, "\t", str(0))
#        output_dataframe = np.vstack([output_dataframe, row])
        print(i, "\t", 0)

# Prep output file:
#output_fname = arg1 + '_tabular'
#np.savetxt(output_fname, output_dataframe, fmt='%s')
sys.stdout.close()
