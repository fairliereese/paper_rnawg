opref=$1
sample=$2
db=$3

d="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/"
sbatch ${d}sbatch_talon.sh $opref $sample $db
qs
