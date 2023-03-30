#!/bin/bash
#SBATCH --job-name=protein_pred
#SBATCH -n 32
#SBATCH -A SEYEDAM_LAB
#SBATCH -o ../processing/%x.o%A
#SBATCH -e ../processing/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=32G
#SBATCH --mail-user=freese@uci.edu

d="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/"

gtf=$1
sample=$2
opref=$(dirname $(dirname $gtf))/protein_pred/$sample
ref=~/mortazavi_lab/data/rnawg/refs/hg38_sirv4_ercc.fa

# convert talon gtf into tama bed
bed=${opref}_tama.bed
tamadir=~/mortazavi_lab/bin/tama/tama_go/format_converter/
python ${tamadir}tama_format_gtf_to_bed12_ensembl.py \
  ${gtf} \
  ${bed}

# convert tama bed into a fasta
module load bedtools2
fasta=${opref}.fa
bedtools getfasta \
  -name \
  -split \
  -s \
  -fi ${ref} \
  -bed ${bed} \
  -fo ${fasta}

# # convert gtf to fasta
# # https://www.biostars.org/p/123166/
# gffdir=~/mortazavi_lab/bin/gffread/
# ${gffdir}./gffread \
#   $gtf \
#   -g $ref \
#   -w $fasta

# call orfs / nmd with tama
tamadir=~/mortazavi_lab/bin/tama/tama_go/orf_nmd_predictions/
orf=${opref}_orf.fa
python ${tamadir}tama_orf_seeker.py \
  -f ${fasta} \
  -o ${orf}

# blast protein sequences against known protein sequences
blastdir=~/mortazavi_lab/bin/ncbi-blast-2.12.0+/bin/
ref_pc=~/mortazavi_lab/data/rnawg/refs/gencode.v29.pc_translations
out=${opref}_blastp.out
${blastdir}./blastp \
  -evalue 1e-10 \
  -num_threads 32 \
  -db ${ref_pc} \
  -query ${orf} > \
  ${out}

# parse blastp output
blastp=$out
parsed_blast=${opref}_blastp_parsed.tsv
python ${tamadir}tama_orf_blastp_parser.py \
  -b ${blastp} \
  -o ${parsed_blast}

# create bed file with cds regions
cds_bed=${opref}_cds.bed
python ${tamadir}tama_cds_regions_bed_add.py \
  -p ${parsed_blast} \
  -a ${bed} \
  -f ${fasta} \
  -o ${cds_bed}

# scan ORFs for protein domains
module load hmmer
out=${opref}_hmmer.out
table=${opref}_hmmer.txt
ref=~/mortazavi_lab/ref/pfam/Pfam-A.hmm
hmmscan \
  -o ${out} \
  --noali \
  --tblout ${table} \
  ${ref} \
  ${orf}
