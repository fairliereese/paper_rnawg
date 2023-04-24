#!/bin/bash

#SBATCH --job-name=suppa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=general,pi_gerstein
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL


module load miniconda
conda activate py3

suppa.py generateEvents -i cerberus.gtf -o localEvents/cerberus.events -f ioe -e SE SS MX RI FL

suppa.py psiPerEvent \
    --ioe-file localEvents/cerberus.events_AF_strict.ioe \
    --expression-file cerberus_filtered_abundance.reformatted.tsv \
    -o psi/cerberus_AF

suppa.py psiPerEvent \
    --ioe-file localEvents/cerberus.events_AL_strict.ioe \
    --expression-file cerberus_filtered_abundance.reformatted.tsv \
    -o psi/cerberus_AL

suppa.py psiPerEvent \
    --ioe-file localEvents/cerberus.events_A3_strict.ioe \
    --expression-file cerberus_filtered_abundance.reformatted.tsv \
    -o psi/cerberus_A3

suppa.py psiPerEvent \
    --ioe-file localEvents/cerberus.events_A5_strict.ioe \
    --expression-file cerberus_filtered_abundance.reformatted.tsv \
    -o psi/cerberus_A5

suppa.py psiPerEvent \
    --ioe-file localEvents/cerberus.events_MX_strict.ioe \
    --expression-file cerberus_filtered_abundance.reformatted.tsv \
    -o psi/cerberus_MX

suppa.py psiPerEvent \
    --ioe-file localEvents/cerberus.events_RI_strict.ioe \
    --expression-file cerberus_filtered_abundance.reformatted.tsv \
    -o psi/cerberus_RI

suppa.py psiPerEvent \
    --ioe-file localEvents/cerberus.events_SE_strict.ioe \
    --expression-file cerberus_filtered_abundance.reformatted.tsv \
    -o psi/cerberus_SE

# For GTF come from GTEx
suppa.py generateEvents \
    -i flair_filter_transcripts.gtf \
    -o localEvents_GTEx/flair.events \
    -f ioe \
    -e SE SS MX RI FL
