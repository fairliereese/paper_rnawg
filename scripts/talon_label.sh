opref=$1
files=$2

# please deal with the damn line endings
dos2unix $files

n=`wc -l $files | cut -d' ' -f1`

d="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/"
sbatch --array=1-${n} ${d}sbatch_talon_label.sh $opref $files
# sbatch --array=29-31 ${d}sbatch_talon_label.sh $opref $files
# sbatch --array=2-9,13,15-18,22-26,30-31,33,37,39-42,45,48-50,52-53,56-66,68-70,73-89 ${d}sbatch_talon_label.sh $opref $files
# sbatch --array=5-6,24-26,50,52,56-57,60-61,65,70,73,75-76,78,80,82-83,85-86 ${d}sbatch_talon_label.sh $opref $files
# sbatch --array=5,25,56,61,75,82-83 ${d}sbatch_talon_label.sh $opref $files
# sbatch --array=83 ${d}sbatch_talon_label.sh $opref $files
