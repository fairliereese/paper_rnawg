opref=$1
ont=$2

d="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/"
sbatch ${d}sbatch_minimap.sh $opref $ont
