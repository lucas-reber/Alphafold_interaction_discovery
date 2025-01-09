# This file includes modifications based on work by Kim et al., 2024, https://doi.org/10.1101/2024.02.19.580970
# from https://github.com/flyark/AFM-LIS/blob/main/alphafold3_lis_contact_v0.15.ipynb, licensed under the MIT License.
#MIT License

#Copyright (c) 2024 Ah-Ram Kim

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


# See the LICENSE file for details.


# The modifications generally change the jupyter notebook to a python file for execution from command line. Only to be used on dimers.
#Usage:
# python LIS_AF3.py -json fold_at1g01550_vs_at2g30680_full_data_0.json
#Output:
# Local Interaction Score_i_avg : 0.000



from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Selection import unfold_entities

import json
import numpy as np
import sys,os
import argparse
import pickle
import itertools
import pandas as pd
import statistics
from collections import Counter
from scipy.optimize import curve_fit

# Modified by Lucas Reber: Allows python file to be executed from command line with specific input
parser = argparse.ArgumentParser(description = '''Calculate LIS from Alphafold3 output ''')

parser.add_argument('-json', nargs=1, type= str, required=True, help ='Input json file.')
parser.add_argument("-cutoff", help="PAE_CUTOFF", nargs='?', type=int, default=12)

def transform_pae_matrix(pae_matrix, pae_cutoff):
    # Initialize the transformed matrix with zeros
    transformed_pae = np.zeros_like(pae_matrix)

    # Apply transformation: pae = 0 -> score = 1, pae = cutoff -> score = 0, above cutoff -> score = 0
    # Linearly scale values between 0 and cutoff to fall between 1 and 0
    within_cutoff = pae_matrix < pae_cutoff
    transformed_pae[within_cutoff] = 1 - (pae_matrix[within_cutoff] / pae_cutoff)
    
    return transformed_pae


def calculate_mean_lis(transformed_pae, subunit_number):
    # Calculate the cumulative sum of protein lengths to get the end indices of the submatrices
    cum_lengths = np.cumsum(subunit_number)

    # Add a zero at the beginning of the cumulative lengths to get the start indices
    start_indices = np.concatenate(([0], cum_lengths[:-1]))

    # Initialize an empty matrix to store the mean LIS
    mean_lis_matrix = np.zeros((len(subunit_number), len(subunit_number)))

    # Iterate over the start and end indices
    for i in range(len(subunit_number)):
        for j in range(len(subunit_number)):
            # Get the start and end indices of the interaction submatrix
            start_i, end_i = start_indices[i], cum_lengths[i]
            start_j, end_j = start_indices[j], cum_lengths[j]

            # Get the interaction submatrix
            submatrix = transformed_pae[start_i:end_i, start_j:end_j]
            

            # Modified by Lucas Reber: Calculate the mean LIS, considering only non-zero values
            non_zero_elements = submatrix[submatrix > 0]
            if non_zero_elements.size > 0:
                mean_lis = non_zero_elements.mean()
            else:
                mean_lis = 0  # or another value, depending on how you want to handle empty submatrices

            # Store the mean LIS in the matrix
            mean_lis_matrix[i, j] = mean_lis

    return mean_lis_matrix


# Modified by Lucas Reber: allows file to be executed from command line and will return the LIS score.
def main():
    pae_cutoff: float = 12.0
    # Your existing code to parse arguments and read the PDB file
    args = parser.parse_args()

    # Initialize the sum of PAE matrices, transformed PAE matrices, mean LIS matrices, and contact maps
    sum_pae_matrix = None
    sum_transformed_pae_matrix = None
    sum_mean_lis_matrix = None
    sum_contact_lia_map = None

    with open(args.json[0], 'r') as file:
        json_data = json.load(file)
    token_chain_ids = json_data['token_chain_ids']
    chain_residue_counts = Counter(token_chain_ids)
    subunit_number = list(chain_residue_counts.values())
    pae_matrix = np.array(json_data['pae'])
    subunit_sizes = subunit_number
    
    # Transform the PAE matrix
    transformed_pae_matrix = transform_pae_matrix(pae_matrix, pae_cutoff)
    transformed_pae_matrix = np.nan_to_num(transformed_pae_matrix)

    # Calculate the mean LIS matrix
    mean_lis_matrix = calculate_mean_lis(transformed_pae_matrix, subunit_sizes)
    mean_lis_matrix = np.nan_to_num(mean_lis_matrix)
    
    value1 = mean_lis_matrix[0, 1]  # Note that indexing in Python starts from 0
    value2 = mean_lis_matrix[1, 0]
    values = [value1, value2]
# Calculate the mean
    if values:
        mean_value = np.mean(values)
    else:
        mean_value = 0
    
    print("Local Interaction Score_i_avg : {:.3f}".format(mean_value))

if __name__ == '__main__':

    main()

