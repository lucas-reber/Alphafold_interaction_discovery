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


#Sripts were modified to allow caculation of LIS for any multimeric complex predicted by Colabfold using Alhpafold-multimer model.
#Usage:
# python LIS_AFM.py -pdb AT1G01550_vs_AT2G30680_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb -json AT1G01550_vs_AT2G30680_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json
#Output:
# Local Interaction Score : 0.000


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

parser = argparse.ArgumentParser(description = '''Calculate LIS of any complex of colabfold output. ''')


parser.add_argument('-pdb', nargs=1, type= str, required=True, help ='Input pdb file.')
parser.add_argument('-json', nargs=1, type= str, required=True, help ='Input pdb file.')
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
    rows, columns = transformed_pae.shape

    # Iterate over the start and end indices
    for i in range(len(subunit_number)):
        for j in range(len(subunit_number)):
            # Get the start and end indices of the interaction submatrix
            start_i, end_i = start_indices[i], cum_lengths[i]
            start_j, end_j = start_indices[j], cum_lengths[j]
            # Get the interaction submatrix
            submatrix = transformed_pae[start_i:end_i, start_j:end_j]
            # Calculate the mean LIS, considering only non-zero values
            

            # Modified by Lucas Reber: Assuming submatrix is a numpy array
            if np.all(submatrix <= 0):
                mean_lis = 0  # Or handle as needed when all elements are <= 0
            else:
                mean_lis = submatrix[submatrix > 0].mean()

            # Store the mean LIS in the matrix
            mean_lis_matrix[i, j] = mean_lis

    return mean_lis_matrix


# Modified by Lucas Reber: Accomodates ColabFold files. allows file to be executed from command line and will return the LIS score
def main():
    pae_cutoff: float = 12.0

    # Initialize the sum of PAE matrices, transformed PAE matrices, mean LIS matrices, and contact maps
    sum_pae_matrix = None
    sum_transformed_pae_matrix = None
    sum_mean_lis_matrix = None
    sum_contact_lia_map = None

    
    ###
    
    args = parser.parse_args()
    pdbp = PDBParser(QUIET=True)
    iopdb = PDBIO()
    structure = pdbp.get_structure('', args.pdb[0])
    
    subunit_sizes = []
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            chain_length = sum(1 for _ in chain.get_residues())
            subunit_sizes.append(chain_length)
    #print(subunit_sizes)
    with open(args.json[0], 'r') as file:
        json_data = json.load(file)
    
    ###

    pae_matrix = np.array(json_data['pae'])
    
    # Transform the PAE matrix
    transformed_pae_matrix = transform_pae_matrix(pae_matrix, pae_cutoff)
    transformed_pae_matrix = np.nan_to_num(transformed_pae_matrix)

    # Calculate the mean LIS matrix
    mean_lis_matrix = calculate_mean_lis(transformed_pae_matrix, subunit_sizes)
    mean_lis_matrix = np.nan_to_num(mean_lis_matrix)
    #print(mean_lis_matrix)
    subunit_number = len(subunit_sizes)
    np.fill_diagonal(mean_lis_matrix, np.nan)

    average_off_diag = np.nanmean(mean_lis_matrix)
    print(f"Local Interaction Score : {average_off_diag:.3f}")


if __name__ == '__main__':

    main()
