### Splitting peak regions into test set (chr2 & chr3) and train set (All chr except chr2 & chr3)

# LR_type="Cerberus"
# LR_type="LAPA"
LR_type=$1

input_dir=data_dir/labeled_beds/all_chr_with_DHS/"$LR_type"/
output_dir=data_dir/labeled_beds/split_chr_with_DHS/"$LR_type"/

train_dir="$output_dir"/train/
test_dir="$output_dir"/test/
test_labels_dir="$output_dir"/test_with_labels/

for inp_bed in "$input_dir"*bed; do
   echo "$inp_bed"
    experiment=$(basename $inp_bed)
    experiment=${experiment##Allchr.}
    experiment=${experiment%%.*}
    cell_line=${experiment%%_ENCSR*}  # GM12878
	echo $experiment
	echo $cell_line
	train_file="$train_dir"/"$experiment".train.bed
	test_file="$test_dir"/"$experiment".test.bed
	test_labels_file="$test_labels_dir"/"$experiment".testlabels.bed
	awk 'BEGIN{FS = OFS = "\t"}{if ($1 != "chr2" && $1 != "chr3") {print $0}}' $inp_bed > $train_file
	awk 'BEGIN{FS = OFS = "\t"}{if ($1 == "chr2" || $1 == "chr3") {NF=NF-2; print}}' $inp_bed > $test_file  # all but last 2 columns = label columns
	awk 'BEGIN{FS = OFS = "\t"}{if ($1 == "chr2" || $1 == "chr3") {print $0}}' $inp_bed > $test_labels_file
done




