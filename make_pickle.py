#Creates pickle file from ColabFold predictions json file, to be used to calculate the pDockQ2 value.
#Usage:
# python make_pickle.py AT1G01550_vs_AT2G30680_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json .
#Output:
# Outputs a output.pkl 


import json
import pickle
import sys
import os

# Check if the input and output directory were provided as command-line arguments
if len(sys.argv) > 2:
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
else:
    print("Usage: python script.py <input_file> <output_directory>")
    sys.exit(1)

# Load data from the specified JSON file
with open(input_file, 'r') as json_file:
    data = json.load(json_file)

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Write the data to a pickle file in the specified output directory
output_file = os.path.join(output_dir, 'output.pkl')
with open(output_file, 'wb') as pickle_file:
    pickle.dump(data, pickle_file)
