#!/bin/bash
#SBATCH --job-name=talon_label
#SBATCH -n 1
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing/%x.o%A_%a
#SBATCH -e processing/%x.e%A_%a
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=16G
#SBATCH --mail-user=freese@uci.edu

# input arguments
opref=$1
samples=$2
genome=~/mortazavi_lab/data/rnawg/refs/hg38_sirv4_ercc.fa

# get the name of the sample
i=$SLURM_ARRAY_TASK_ID
p=`head -${i} $samples | tail -1`

sam=${opref}${p}.bam
echo $sam
talon_label_reads \
    --f $sam \
    --g $genome \
    --tmpDir ${opref}${p}_tmp/ \
    --ar 20  \
    --deleteTmp  \
    --o ${opref}${p}
