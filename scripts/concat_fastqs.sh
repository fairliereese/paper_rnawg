files=$1
opref=$2

# please deal with the damn line endings
dos2unix $files

d="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/"
sbatch ${d}sbatch_concat_fastqs.sh $files $opref
