#!/bin/bash
#SBATCH --job-name=talon_annot
#SBATCH -n 16
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing/%x.o%A
#SBATCH -e processing/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=128G
#SBATCH --mail-user=freese@uci.edu

# usage
# sbatch sbatch_talon_read_annot.sh <oprefix>

opref=$1
build=hg38
db=${opref}.db

talon_fetch_reads \
    --db ${db} \
    --build $build \
    --o ${opref}
