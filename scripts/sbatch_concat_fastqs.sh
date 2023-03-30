#!/bin/bash
#SBATCH --job-name=concat_fastqs
#SBATCH -n 32
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing/%x.o%A
#SBATCH -e processing/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=4:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=64G
#SBATCH --mail-user=freese@uci.edu

# input tab-separated file with PB ID \t file location at GHTF \t md5sum
ifile=$1
opref=$2

ofile=${opref}.fastq
touch ${ofile}
rm $ofile
touch ${ofile}

while read sample
do
  # extract PBID
  pb_id=`echo $sample | cut -f1 -d' '`
  smrt_cell=`echo $sample | cut -f4 -d' '`

  # make directories
  pb_dir=~/pacbio/$pb_id/

  # get flnc post-Refine reads for each data directory
  dir=${pb_dir}Refine/${smrt_cell}01/
  files=($dir/flnc.fastq)
  fastq=${files[0]}
  cat ${fastq} >> ${ofile}
done < ${ifile}
