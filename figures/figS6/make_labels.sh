### Assigns binary labels to peak regions reported by LR assays and methods (Cerb/LAPA)
# based on the peak presence in other assays (i.e RAMPAGE & CAGE)

# LR_type="Cerberus"
# LR_type="LAPA"
LR_type=$1

input_dir=data_dir/processed_beds/"$LR_type"/
output_dir=data_dir/labeled_beds/all_chr/"$LR_type"/
inp_assay_dir=data_dir/all_chr/

for inp_bed in "$input_dir"*.bed; do
		echo $inp_bed
    experiment=$(basename $inp_bed)
    experiment=${experiment%%.bed}
	cell_line=${experiment%%_ENCSR*}  # GM12878
    assay1="CAGE"
    assay2="RAMPAGE"
    echo $experiment
    echo $cell_line
	echo "------"
    for other_assay in $assay1 $assay2; do
        other_assay_file="$inp_assay_dir"/Allchr."$other_assay"-"$cell_line".txt
        tmp_other_assay_intersect="$output_dir"/"$experiment"."$other_assay".intersect.tmp
        tmp_other_assay_notintersect="$output_dir"/"$experiment"."$other_assay".notintersect.tmp
        bedtools intersect -wa -a <(sort -k1,1 -k2,2n $inp_bed) -b "$other_assay_file" | sort -k1,1 -k2,2n | uniq > $tmp_other_assay_intersect
        bedtools intersect -wa -v -a <(sort -k1,1 -k2,2n $inp_bed) -b $other_assay_file | sort -k1,1 -k2,2n | uniq > $tmp_other_assay_notintersect
        awk -v assay="$other_assay" 'BEGIN{FS = OFS ="\t"}{print $0, assay}' \
                $tmp_other_assay_intersect > $tmp_other_assay_intersect.tmp
        awk -v assay="---" 'BEGIN{FS = OFS ="\t"}{print $0, assay}' \
                $tmp_other_assay_notintersect > $tmp_other_assay_notintersect.tmp
        tmp_other_assay_file="$output_dir"/"$experiment"."$other_assay".tmp
        cat $tmp_other_assay_intersect.tmp $tmp_other_assay_notintersect.tmp | sort -k1,1 -k2,2n > $tmp_other_assay_file
    done
    tmp_other_assay1="$output_dir"/"$experiment"."$assay1".tmp
    tmp_other_assay2="$output_dir"/"$experiment"."$assay2".tmp
    output_labeled_file="$output_dir"/Allchr."$experiment".labeled.bed
    paste -d ',' "$tmp_other_assay1" <(cut -f 6 $tmp_other_assay2) | awk 'BEGIN{FS = OFS = "\t"} {if ($NF=="---,---") {print $0, 0} else {print $0, 1} }' > $output_labeled_file
    rm "$output_dir"/*.tmp
done



