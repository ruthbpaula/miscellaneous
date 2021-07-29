# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 18:56:06 2021

@author: Ruth De Paula and Luciano Branco
"""

import pandas as pd
import numpy as np
import math
import numpy.matlib

# # arg1 = sys.argv[1]
arg1 = 'all_g4_intron_single_coord' # input file with split coordinates

file_3utr_genes = pd.read_csv(arg1, delimiter = '\t', names = ['chro', 'coord']) # ficticious header to access columns with (.column_name)
file_3utr_genes_numpy = file_3utr_genes.to_numpy()

# we need a first comparison to start the loop
key_for_comparison = file_3utr_genes_numpy[0]
possible_endpoint = key_for_comparison[1]
index = 0

# remove first element since it is already saved
file_3utr_genes_numpy = np.delete(file_3utr_genes_numpy, 0, 0)
output_dataframe = np.empty([0,3])

for current_line in file_3utr_genes_numpy: # current_line = file_3utr_genes_numpy[6] 

    index += 1
    
    if (current_line[0] == key_for_comparison[0]):
        if (current_line[1] == (key_for_comparison[1] + index)):
            # possible match. Verify next line
            possible_endpoint = current_line[1]
        else:
            # match found! save the output row and edit the key
            row = np.append(key_for_comparison, possible_endpoint)
            output_dataframe = np.vstack([output_dataframe, row])
            key_for_comparison = current_line
            index = 0
        
    else: # chr name has changed. match
    
        row = np.append(key_for_comparison, possible_endpoint)
        output_dataframe = np.vstack([output_dataframe, row])
        key_for_comparison = current_line
        index = 0

row = np.append(key_for_comparison, possible_endpoint)
output_dataframe = np.vstack([output_dataframe, row])

# Prep output file:
output_fname = arg1 + '_tabular'
np.savetxt(output_fname, output_dataframe, fmt='%s')
