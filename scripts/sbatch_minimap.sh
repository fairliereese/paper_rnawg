#!/bin/bash
#SBATCH --job-name=minimap
#SBATCH -n 32
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing/%x.o%A
#SBATCH -e processing/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=64G
#SBATCH --mail-user=freese@uci.edu

opref=$1
ont=$2

fastq=${opref}_demux.fastq
ref=~/mortazavi_lab/data/rnawg/refs/hg38_sirv4_ercc.fa
sj_ref=~/mortazavi_lab/ref/gencode.vM21/gencode.vM21.bed
sam=${opref}_mapped.sam
log=${opref}_minimap.log


module load minimap2

if [ -z "$ont" ]
  then
    # ont settings - alternative
    minimap2 \
      -t 32 \
      -ax splice \
      -uf \
      --MD \
      --secondary=no \
      --junc-bed $sj_ref \
      $ref \
      $fastq \
      > $sam
      2> $log
  else
    # pacbio settings - default
    minimap2  \
        -t 32 \
        -ax splice:hq \
        -uf \
        --MD \
        $ref \
        $fastq \
        > $sam \
        2> $log
fi
