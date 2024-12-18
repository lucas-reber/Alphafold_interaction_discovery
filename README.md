# Alphafold_interaction_discovery

This repository provides details on using AlphaFold for protein-protein interaction discovery and was used for analysis in the paper XXX (unpublished).
We used AlphaFold-Multimer (as implemented in [ColabFold](https://doi.org/10.1038/s41592-022-01488-1)) and AlphaFold3 (using the [AlphaFold Server](https://golgi.sandbox.google.com/about)) for predictions of interactions.
We used different metrics to evaluate predictions: model confidence (pTM + ipTM), ipTM, pDockQ, pDockQ2 and LIS (Bryant et al., 2022; Evans et al., 2022; Jumper et al., 2021; A.-R. Kim et al., 2024; Zhu et al., 2023). Details on how to extract these metrics from AlphaFold-Multimer and AlphaFold3 outputs can be found below. 

# Running predictions
[Details on how to install ColabFold](https://github.com/YoshitakaMo/localcolabfold)  

We ran predictions through ColabFold ([Mirdita et al., 2022](https://doi.org/10.1038/s41592-022-01488-1)) in two steps. First the input features were generated using the `colabfold_search` and then the predictions were run on the input features using `colabfold_batch`. This approach allows to perform the complete analysis locally without connection to the ColabFold server and therefore is not limited by the ColabFold server. However, this requires setting up local databases for genetic and template search. 
For information on the databases please visit the [ColabFold repository](https://github.com/sokrypton/ColabFold/tree/main). This also requires the PDB templates similar to the original AlphaFold2 implementation.  

A detailed `config.json` with the parameters we used for ColabFold predictions can be found in `example_files`.
## Colabfold_search - Input feature generation
This step requires a lot of RAM and CPUs. We used 16 CPUs (Intel(R) Xeon(R) Gold 6448H) and 280GB of RAM. With these specifications it took 6-12 h for 5000 sequences.
The input file is a csv file formatted like the `EXAMPLE_INPUT.csv` file.
For multimer predictions the multiple protein sequences are all pasted in one line but separated by a colon ":".

Further requirements:  
requires conda installation (e.g. Miniconda3/22.11.1-1)  
requires mmseqs2 (e.g. MMseqs2/15-6f452)  
requires database to be setup (https://colabfold.mmseqs.com)  


```
# activate the colabfold environment
eval "$(conda shell.bash hook)"
conda activate /path/to/localcolabfold/colabfold-conda

colabfold_search \
  --use-env 1 \
  --use-templates 1 \ #will generate .m8 template features
  --db-load-mode 2 \
  --mmseqs /path/to/bin/mmseqs \
  --db2 pdb100_230517 \
  --threads 16 \
  EXAMPLE_INPUT.csv \ #input csv file
  /path/to/database/mmseqs/uniref30_2302 \ 
  msas #output directory
```
## Colabfold_batch - Structure prediction

This step requires a GPU. We used 80 GB A100 or H100 nVidia GPUs. A L40S 48GB GPU was also successfully used. In addtion we used 2 CPU cores (Intel(R) Xeon(R) Gold 6326 CPU @ 2.90GHz) and 70-150GB of RAM depending on protein length. Proteins of over 4000 amino acids were the approximate limit.
A single GPU can do approximately 50 prediction per day depending on protein length. We frequently used 30-60 GPUs in parallel allowing us to do around 1500-3000 predictions in a day. Note that this script should be run on an individual GPU, but multiple instances may be run at the same time.

Further requirements:  
requires conda installation (e.g. Miniconda3/22.11.1-1)  
requires CUDA driver (e.g. CUDA/12.2.0)
requires the AlphaFold2 templates (see https://github.com/google-deepmind/alphafold)  

General setup:  
```
eval "$(conda shell.bash hook)"

# activate the colabfold environment
conda activate /path/to/localcolabfold/colabfold-conda

# Define the local PDB path
LOCALPDBPATH="/path/to/v2.3.2/pdb_mmcif/mmcif_files"

# Define the random seed
RANDOMSEED=0
```
Run prediction on a single input file:  
```
#specify input files
INPUTFILE="input_file.a3m"
PDBHITFILE="input_file_pdb100_230517.m8"

#check if .m8 is empty 
if [ -s "${PDBHITFILE}" ]; then
    
    colabfold_batch \
      --amber \ 
      --templates \
      --use-gpu-relax \
      --pdb-hit-file ${PDBHITFILE} \
      --local-pdb-path ${LOCALPDBPATH} \
      --random-seed ${RANDOMSEED} \
      ${INPUTFILE} \
      output_dir
   
  else
    
    colabfold_batch \
      --amber \
      --use-gpu-relax \
      --random-seed ${RANDOMSEED} \
      ${INPUTFILE} \
      output_dir

fi
```
Run prediction on multiple files in a directory:  
```
# Loop over each .a3m file in the msas directory
for INPUTFILE in "${msas_folder}"/*.a3m; do

  # Extract the base name from the INPUTFILE name
  NAME=$(basename "$INPUTFILE" .a3m)

  # Define the PDB hit file
  PDBHITFILE="${msas_folder}/${NAME}_pdb100_230517.m8"
  
  # Define prediction directory for the output
  pred_dir="pred_${NAME}"

  #check if the .m8 file has any hits.
  if [ -s "${PDBHITFILE}" ]; then
    # Run the colabfold_batch command
    colabfold_batch \
      --amber \
      --templates \
      --use-gpu-relax \
      --pdb-hit-file ${PDBHITFILE} \
      --local-pdb-path ${LOCALPDBPATH} \
      --random-seed ${RANDOMSEED} \
      ${INPUTFILE} \
      ${pred_dir}

  else
    # Run the colabfold_batch command
    colabfold_batch \
      --amber \
      --use-gpu-relax \
      --random-seed ${RANDOMSEED} \
      ${INPUTFILE} \
      ${pred_dir}
  fi

done
```
# Extracting confidence metrics
By default AlphaFold will produce five predictions which are based on five different models. Using all the five outputs and averaging metrics across the outputs can improve the performance of the metrics for the selection of interaction candidates.
## ipTM and pTM

The interface predicted TM score (ipTM) and predicted TM score (pTM) are directly provided by AlphaFold predictions. They can be found in the json file (e.g. in `example_files` for Colabfold `AT1G01550_vs_AT2G30680_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json` and for AlphaFold3 `fold_at1g01550_vs_at2g30680_summary_confidences_0.json`).  


We extracted the scores using the following bash code:

```
  file="AT1G01550_vs_AT2G30680_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json"
  ptm_score=$(jq '.ptm' "$file")
  iptm_score=$(jq '.iptm' "$file")
  echo "$ptm_score,$iptm_score"
```

To calculate the model confidence score use this formula: Model confidence = ipTM * 0.8 + pTM * 0.2.

## LIS
The Local Interaction Score was developed by [A.-R. Kim et al., 2024](https://doi.org/10.1101/2024.02.19.580970).
We obtained scripts from the [AFM-LIS repository](https://github.com/flyark/AFM-LIS) and modified them where needed.
Use `LIS-AFM.py` script to extract LIS score from AlphaFold-Multimer predictions performed through Colabfold.  
Usage:
```
python LIS_AFM.py -pdb AT1G01550_vs_AT2G30680_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb -json AT1G01550_vs_AT2G30680_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json
```
Output:
```
Local Interaction Score : 0.000
```

Use `LIS-AF3.py` script to extract LIS score from AlphaFold2 predictions performed through the AlphaFold Server.
Usage:
```
python LIS_AF3.py -json fold_at1g01550_vs_at2g30680_full_data_0.json
```
Output:
```
Local Interaction Score_i_avg : 0.000
```




## pDockQ
Only applicable for AlphaFold-Multimer/ColabFold predictions.
We used [pDockQ scripts](https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/) by [Bryant et al., 2022](https://doi.org/10.1038/s41467-022-28865-w)

Usage: 
```
python pdockq.py --pdbfile AT1G01550_vs_AT2G30680_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb
```
Output: 
```
pDockQ = 0.058 for AT1G01550_vs_AT2G30680_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb
This corresponds to a PPV of at least 0.63555449
```
## pDockQ2
Only applicable for AlphaFold-Multimer/ColabFold predictions.
We used [pDockQ2 scripts](https://gitlab.com/ElofssonLab/afm-benchmark/-/blob/main/src/pdockq2.py) by [Zhu et al., 2023](https://doi.org/10.1093/bioinformatics/btad424)  
This scripts takes a pickle file as input which was part of the output of the original AlphaFold2 implementation but not part of Colabfold output.  
We provide a `make_pickle.py` script to create pickle files from ColabFold json files.
Usage:
```
python make_pickle.py AT1G01550_vs_AT2G30680_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json .
```
The output is a `output.pkl` file and can be used in the next step to extract pDockQ2 scores.

Usage: 
```
python pdockq2.py -pkl output.pkl -pdb AT1G01550_vs_AT2G30680_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb -dist 10
```
Output: 
```
pDockQ_i is:
A 0.010403525668536952
B 0.008154758661146998
```
We used the average of the per chain pDockQ_A and pDockQ_B to calculate the pDockQ2.
# Example calculations on metrics
The `example_calculations.R` script allows to perform perform basic calculations such as deteriming the best and average metrics across AlphaFold-Multimer and AlphaFold3 outputs, and calculating model confidence (0.2*pTM + 0.8*ipTM) or pDockQ2. It also shows how to calculate the Mean AFM score and the AFM & AF3 score from those metrics. The script can be run with `example_data_1.csv` and `example_data_2.csv`. This script was used to generate Dataset S1, Dataset S2, Dataset S3, Dataset S4 and Dataset S5 from our publication.


# ROC analysis
The `ROC_analysis.R` script is provided and allows to perform the Receiver Operating Characteristic (ROC) curves and Area Under the Curve (AUC) statistics seen in our paper in SI Appendix Fig. S5D (to be published).
Example data is deposited in `example_files/example_data_3.csv` and the full analysis can be run on Dataset S1 from the publication.



# References
[Run ColabFold on your local computer](https://github.com/YoshitakaMo/localcolabfold) by Yoshitaka Moriwaki  
[Bryant et al., 2022](https://doi.org/10.1038/s41467-022-28865-w)  
[Zhu et al., 2023](https://doi.org/10.1093/bioinformatics/btad424)  
[A.-R. Kim et al., 2024](https://doi.org/10.1101/2024.02.19.580970)  
[Jumper et al., 2021](https://doi.org/10.1038/s41586-021-03819-2)  
[Evans et al., 2022](https://doi.org/10.1101/2021.10.04.463034)  
[Mirdita et al., 2022](https://doi.org/10.1038/s41592-022-01488-1)  
[Abramson et al., 2024](https://doi.org/10.1038/s41586-024-07487-w)

