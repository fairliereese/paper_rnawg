### Calculating DNase-seq (DHS) signal values over each peak region
### Combining peak location, peak id, assay expression, DHS and overlap labels into one file

# LR_type="Cerberus"
# LR_type="LAPA"
LR_type=$1

input_dir=data_dir/labeled_beds/all_chr/"$LR_type"/
output_dir=data_dir/labeled_beds/all_chr_with_DHS/"$LR_type"/
bw_dir=data_dir/bigWigs/

for inp_bed in "$input_dir"*bed; do
	experiment=$(basename $inp_bed)
    experiment=${experiment##Allchr.}
    experiment=${experiment%%.bed}
    cell_line=${experiment%%_ENCSR*}  # GM12878
	echo $experiment
	echo $cell_line
	echo "-------"
	DHS_bw="$bw_dir"/"$cell_line".DHS.bigWig
	bw_output="$output_dir"/"$experiment".DHS.tab
	./bigWigAverageOverBed $DHS_bw <(awk 'BEGIN{FS=OFS="\t"}{print $1, $2, $3, $4}' $inp_bed) $bw_output
	
	experiment2=${experiment%%.labeled}
	out_bed="$output_dir"$experiment2.DHS.bed
	# output format: chr	start	end		peak_id	TPM	DHS	full_label(CAGE,---)	label(0/1)
	join -1 4 -2 1 -o 1.1,1.2,1.3,1.4,1.5,2.6,1.6,1.7 <(sort -k4,4 $inp_bed) <(sort -k1,1 $bw_output) | sed "s# #\t#g" | \
			sort -k1,1 -k2,2n > $out_bed
	rm "$output_dir/"*DHS.tab
done



