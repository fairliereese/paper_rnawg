## Input files
1. `swan_metadata.txt`
2. `cerberus_psi.tsv`
3. `suppa/psi/cerberus_AF.psi` and `suppa/psi/cerberus_AL.psi` from the SUPPA output (see S7A).
4. `lr_human_library_data_summary.tsv`

## Compare AS genes identified by Cerberus and/or SUPPA
Execute `cerberus_suppa.ipynb` to obtain AS genes identified by Cerberus and/or SUPPA with default thresholds (0.25-0.75).

## Visualization
Run `comp_cerberus_suppa_with_different_thres.R` to generate `cerberus_suppa_psi_0.25_0.75.pdf`.