# coding=utf-8
##intersecting files

import pandas as pd
import numpy as np
import math
import numpy.matlib


# # arg1 = sys.argv[1]
arg1 = 'gross-dels_toIntersect_ok' # input file with HS coordinates

# # arg2 = sys.argv[2]
arg2 = 'chrs_endsize.txt' # input file with chr sizes


file_hyper_sens = pd.read_csv(arg1, delimiter = '\t', names = ['chr', 'start', 'end']) # ficticious header to access columns with (.column_name)
file_size_chrs = pd.read_csv(arg2, delimiter = '\t', names = ['chr', 'size_chr'])
file_size_chrs_numpy = file_size_chrs.to_numpy()

bin_size = 500

# my_chro[0] is chro name ; my_chro[1] is chro maxsize
# my_chro =file_size_chrs_numpy[0]
for my_chro in file_size_chrs_numpy:
    chro_name   = my_chro[0]
    chro_length = my_chro[1]
    print('Searching ' + str(chro_name) + ', length ' + str(chro_length) + '...')
    
    num_bins = math.floor(chro_length / bin_size)
    
    nucleotids_occurs = np.zeros(num_bins, dtype=int) # Create zero positions
    
    # -- Search hsens file and find matching occurences: --
    
    # Get all chros from sens file, matching my_chro:
        # DropNA remove rows not matching 'my_chro[0]'
        # to_numpy() for performance
    chro_sens = file_hyper_sens.where(file_hyper_sens.chr == chro_name).dropna().to_numpy()
    
    # hsens =chro_sens[0]
    for hsens in chro_sens:
        hs_start = hsens[1]
        hs_end   = hsens[2]
        bin_start = math.ceil(hs_start / bin_size) -1
        bin_end   = math.ceil(hs_end / bin_size) -1
        
        nucl_occur = hs_end - hs_start
        
        if (bin_start >= num_bins or bin_end >= num_bins):
            error_msg = '-> '+ str(chro_name) + ' - Skipped trying to access ' + str(bin_start) + ':' + str(bin_end) + ' of chr size ' + str(num_bins) + '. Please verify input files.'
            print(error_msg)
            break # Avoid error out of bounds for 'nucleotids_occurs' vector
        
        if (bin_end == bin_start):
            nucleotids_occurs[bin_end] += nucl_occur
        else:
            add_to_first_bin = bin_size * (bin_start+1) - hs_start +1
            nucleotids_occurs[bin_start] += add_to_first_bin
            
            nucl_occur -= add_to_first_bin # Subtract whats already accounted for
            
            nucleotids_occurs[bin_end] += nucl_occur - (bin_end - bin_start -1)*500
            
            if (bin_end - bin_start > 1):
                nucleotids_occurs[bin_start+1 : bin_end] += 500
    
    # Prep output dataframe
    out_col1 = np.squeeze( numpy.matlib.repmat(chro_name, 1, num_bins) )
    out_col2 = np.arange(1, chro_length-bin_size+1, bin_size, dtype=int)
    out_col3 = np.arange(bin_size+1, chro_length+1, bin_size, dtype=int)
    
    chro_out = {'chr':out_col1 , 'bin_sta': out_col2, 'bin_end': out_col3, 'occrs':nucleotids_occurs}
    chro_out_df = pd.DataFrame(chro_out)
    
    # Prep output file:
    out_fname = chro_name + '_occur.txt'
    chro_out_df.to_csv(out_fname, header=False, index=False, sep='\t', mode='w')
    

del file_hyper_sens
del file_size_chrs
# output_file.close()

# # my_chro[0] is chro name ; my_chro[1] is chro maxsize
# for my_chro in file_size_chrs_numpy:
#     chro_length = my_chro[1]
#     for i in range(1, chro_length, bin_size):
#         # Searching from i to bin_size-1 ...
        
#         # Get all chros from sens file, matching my_chro:
#             # DropNA remove rows not matching 'my_chro[0]'
#             # to_numpy() for performance
#         chro_sens = file_hyper_sens.where(file_hyper_sens.chr == my_chro[0]).dropna().to_numpy()
        
        
#         print('From: ', i)
#         print('To: ',i+bin_size-1)
        
#         if (i > 3):
#             break
    
    
#     break


