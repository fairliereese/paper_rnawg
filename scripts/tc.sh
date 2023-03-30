opref=$1

d="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/"
sbatch ${d}sbatch_tc.sh $opref
