#!/bin/bash
#SBATCH --job-name=talon
#SBATCH -n 16
#SBATCH -A SEYEDAM_LAB
#SBATCH -o processing/%x.o%A
#SBATCH -e processing/%x.e%A
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END
#SBATCH --mem=256G
#SBATCH --mail-user=freese@uci.edu

# usage
# sbatch sbatch_talon_bulk.sh -d <database> <config> <oprefix>

# get -d (database) args
while getopts d: flag; do
  case "$flag" in
    d) db=${OPTARG} ;;
  esac
done

# shift flag args away and get opref / sample
shift $((OPTIND - 1))

config=$1
opref=$2

gtf=~/mortazavi_lab/data/rnawg/refs/gencode_v29_sirv4_ercc.gtf
build=hg38

if [ -z "$db" ]
  then
    echo "No database given. Will make new database."

    talon_initialize_database \
        --f ${gtf} \
        --g ${build} \
        --a gencode_v29 \
        --l 0 \
        --idprefix ENCODEH \
        --5p 500 \
        --3p 300 \
        --o ${opref}

      db=${opref}.db

  else
    echo "Adding reads to database ${db}"
fi

talon \
      --f ${config} \
      --db ${db} \
      --build $build \
      -t 32 \
      --o ${opref}
