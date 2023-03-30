#!/bin/bash
#SBATCH --job-name=demux
#SBATCH -n 64
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing/%x.o%A
#SBATCH -e processing/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=48:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=256G
#SBATCH --mail-user=freese@uci.edu

opref=$1
fastq=${opref}.fastq
python ~/mortazavi_lab/bin/LR-splitpipe/LR-splitpipe/demultiplex.py \
    -f ${fastq} \
    -o ${opref} \
    -t 32 \
    -rc 0
