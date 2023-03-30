#!/bin/bash
#SBATCH --job-name=talon_ab_gtf
#SBATCH -n 1
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing/%x.o%A
#SBATCH -e processing/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=128G
#SBATCH --mail-user=freese@uci.edu

set -x
set -e

datasets=false

# get -d (database) args
while getopts d: flag; do
  case "$flag" in
    d) datasets=${OPTARG} ;;
  esac
done

# shift flag args away and get opref / sample
shift $((OPTIND - 1))

db=$1
opref=$2

annot=gencode_v29
build=hg38

if [ ! "$datasets" ]
then
  # unfiltered talon abundance
  # used for gene level quantification
  talon_abundance \
      --db $db \
      -a ${annot} \
      -b ${build} \
      -d ${datasets} \
      --o ${opref}

  # filter transcripts
  talon_filter_transcripts \
      --db ${db} \
      -a ${annot} \
      --datasets ${datasets} \
      --maxFracA 0.5 \
      --minCount 5 \
      --minDatasets 2 \
      --o ${opref}_pass_list.csv

  # filtered talon abundance
  # used for transcript level quantification
  talon_abundance \
      --db ${db} \
      -a ${annot} \
      -b ${build} \
      -d ${datasets} \
      --whitelist ${opref}_pass_list.csv \
      --o ${opref}

  # filtered GTF
  talon_create_GTF \
      --db ${db} \
      -a ${annot} \
      -b ${build} \
      -d ${datasets} \
      --whitelist ${opref}_pass_list.csv \
      --o ${opref}

  # get pass list with only known, NIC, NNC
  d=~/mortazavi_lab/data/rnawg/scripts/
  ab=${opref}_talon_abundance_filtered.tsv
  python ${d}get_known_nic_nnc_pass_list.py \
    ${ab} \
    ${opref}_pass_list.csv \
    ${opref}

  # filtered GTF with only known, NIC, NNC
  talon_create_GTF \
      --db ${db} \
      -a ${annot} \
      -b ${build} \
      -d ${datasets} \
      --whitelist ${opref}_known_nic_nnc_pass_list.csv \
      --o ${opref}_known_nic_nnc

  else
    # unfiltered talon abundance
    # used for gene level quantification
    talon_abundance \
        --db $db \
        -a ${annot} \
        -b ${build} \
        --o ${opref}

    # filter transcripts
    talon_filter_transcripts \
        --db ${db} \
        -a ${annot} \
        --maxFracA 0.5 \
        --minCount 5 \
        --minDatasets 2 \
        --o ${opref}_pass_list.csv

    # filtered talon abundance
    # used for transcript level quantification
    talon_abundance \
        --db ${db} \
        -a ${annot} \
        -b ${build} \
        --whitelist ${opref}_pass_list.csv \
        --o ${opref}

    # filtered GTF
    talon_create_GTF \
        --db ${db} \
        -a ${annot} \
        -b ${build} \
        --whitelist ${opref}_pass_list.csv \
        --o ${opref}

    # get pass list with only known, NIC, NNC
    d=~/mortazavi_lab/data/rnawg/scripts/
    ab=${opref}_talon_abundance_filtered.tsv
    python ${d}get_known_nic_nnc_pass_list.py \
      ${ab} \
      ${opref}_pass_list.csv \
      ${opref}

    # filtered GTF with only known, NIC, NNC
    talon_create_GTF \
        --db ${db} \
        -a ${annot} \
        -b ${build} \
        --whitelist ${opref}_known_nic_nnc_pass_list.csv \
        --o ${opref}_known_nic_nnc
  fi
