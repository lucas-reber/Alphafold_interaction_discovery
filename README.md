# Alphafold_interaction_discovery

The code was used for analysis in the paper XXX.
We used AlphaFold-Multimer (through Colabfold) and AlphaFold3 to screen for protein-protein interactions.
To validate this approach we evaluated different confidence metrics on a arabidopsis protein-protein interaction reference set.
Model confidence (pTM + ipTM), ipTM, pDockQ, pDockQ2 and LIS (Bryant et al., 2022; Evans et al., 2022; Jumper et al., 2021; A.-R. Kim et al., 2024; Zhu et al., 2023). This repository includes details on how to obtain the metrics from Colabfold and AlphaFold3 outputs. By default AlphaFold will produce five predictions which are based on five different models. Using all the five outputs and averaging scores across the outputs can improve the performance of the scores.

# Running predictions

Colabfold

# Extracting confidence metrics

## ipTM and pTM

The interface predicted TM score (ipTM) and predicted TM score (pTM) are directly provided by AlphaFold predictions. They can be found in the json file (e.g. in `example_files` for Colabfold `AT1G01550_vs_AT2G30680_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json` and for AlphaFold3 `fold_at1g01550_vs_at2g30680_summary_confidences_0.json`). 
We extracted the scores using bash using following code:
```
  file="AT1G01550_vs_AT2G30680_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json"
  ptm_score=$(jq '.ptm' "$file")
  iptm_score=$(jq '.iptm' "$file")
  echo "$ptm_score,$iptm_score"
```

To calculate the model confidence score use this formula: Model confidence = ipTM * 0.8 + pTM * 0.2.
