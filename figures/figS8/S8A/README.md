## Input files
1. `suppa/cerberus.gtf`
2. `suppa/flair_filter_transcripts.gtf`
3. `suppa/cerberus_filtered_abundance.tsv`
4. `cerberus_transcript_novelty.tsv`

## Data Pre-processing 
1. Reformat the GTF file: `sed -i 's/,/_/g' cerberus.gtf`.
2. Filter the abundance matrix based on transcripts in the novelty file: `Rscript filter_abundance_mat_by_novelty.R`.
3. Reformat the abundance matrix as the input of the SUPPA: `cut -f4,12- cerberus_filtered_abundance.novelty.tsv > tmp.tsv ; sed -i 's/,/_/g' tmp.tsv; mv tmp.tsv cerberus_filtered_abundance.reformatted.tsv`.
4. Remove the name of first column manually, as required by SUPPA.

## Run SUPPA
To generate local events based on `cerberus.gtf` and ` flair_filter_transcripts.gtf`, and calculate PSI values for local events, please run `sbatch suppa/suppa.slurm`.
1. `suppa/localEvents`: local events generated based on `cerberus.gtf`.
2. `suppa/localEvents_GTEx`: local events generated based on `flair_filter_transcripts.gtf`.
3. `suppa/psi`: PSI values calculated for local events generated based on `cerberus.gtf` using expression matrix `cerberus_filtered_abundance.reformatted.tsv`.

## Generate Fig. S7A
Run `Rscript prop_of_novel_txs_in_localEvents.R` to produce the raw Fig. S7A (`proportion.pdf`), displaying the proportion of observed known and novel transcripts identified by Cerberus, based on novel events at the TSS, TES, or EC. Note, the total number of local events is shown in the file `suppa/tot_number_of_localEvents.tsv`.