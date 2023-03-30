#!/bin/bash
#SBATCH --job-name=add_bam_tag
#SBATCH -n 1
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing/%x.o%A
#SBATCH -e processing/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=16G
#SBATCH --mail-user=freese@uci.edu

opref=$1

sam=${opref}_clean.sam

python ~/mortazavi_lab/bin/LR-splitpipe/LR-splitpipe/add_bam_tag.py \
  -s $sam \
  --merge_primers \
  -o ${opref}
