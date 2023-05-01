## Machine learning models predict the support for long-read TSS peaks by other TSS-annotating assays and in a cross-cell type manner

# After running this script, Supp. Fig. 6 panels should be created in the following directories:
# SFig 6a: plots/encode_rampage_cage.pdf
# SFig 6b: plots/overlap.pdf
# SFig 6c: plots/test_train_set_peaks.pdf
# SFig 6d: training/plots/Cerberus_AIC.png  training/plots/LAPA_AIC.png
# SFig 6e: training/plots/individual_exper_same_cell_line/Cerberus/GM12878_ENCSR962BVU_AUROC.pdf training/plots/individual_exper_same_cell_line/LAPA/GM12878_ENCSR962BVU_AUROC.pdf
# SFig 6f: plots/AUROC_diff_cell_line_boxplot.pdf

#########################################################################

# Downloading DHS-seq bigWig files from ENCODE for GM12878 and K562
mkdir -p data_dir/bigWigs
curl -J -O -L "https://www.encodeproject.org/files/ENCFF428XFI/@@download/ENCFF428XFI.bigWig"
curl -J -O -L "https://www.encodeproject.org/files/ENCFF414OGC/@@download/ENCFF414OGC.bigWig"
mv ENCFF428XFI.bigWig data_dir/bigWigs/GM12878.DHS.bigWig
mv ENCFF414OGC.bigWig data_dir/bigWigs/K562.DHS.bigWig
echo "############################# .bigwig download is done"

# Cerberus replicates will be combined into a file in this directory
mkdir -p data_dir/raw_beds/Cerberus 
Rscript combine_cerberus_replicates.R
echo "############################# Cerberus replicates are combined"

# We have raw_beds/LAPA and raw_beds/Cerberus bed files for the LR experiments. Sorting peaks and filtering out any chromosomes other than chr1-22,XY:
mkdir -p data_dir/processed_beds/Cerberus data_dir/processed_beds/LAPA
Rscript preproc_files.R Cerberus
Rscript preproc_files.R LAPA
echo "############################# Sorting peaks and preprocessing are done"

# We have the files in processed_beds/*/*bed. Now creating binary labels based on LR peak presence in RAMPAGE/CAGE assays
mkdir -p data_dir/labeled_beds/all_chr/Cerberus data_dir/labeled_beds/all_chr/LAPA
./make_labels.sh Cerberus
./make_labels.sh LAPA
echo "############################# Binary labels created based on peak overlap with RAMPAGE and CAGE assays."

# Finding DHS values over peaks
mkdir -p data_dir/labeled_beds/all_chr_with_DHS/Cerberus data_dir/labeled_beds/all_chr_with_DHS/LAPA
./bigwigavg_DHS.sh Cerberus
./bigwigavg_DHS.sh LAPA
echo "############################# DHS signal over peaks is calculated."

# Creating test and train set (chr2/3 vs. other chrs)
for LR in Cerberus LAPA; do
	mkdir -p data_dir/labeled_beds/split_chr_with_DHS/"$LR"
	mkdir -p data_dir/labeled_beds/split_chr_with_DHS/"$LR"/train
	mkdir -p data_dir/labeled_beds/split_chr_with_DHS/"$LR"/test
	mkdir -p data_dir/labeled_beds/split_chr_with_DHS/"$LR"/test_with_labels
	./split_train_test.sh "$LR"
done
echo "############################# train/test sets are created. Changing the current directory to 'training' ..."

# Preprocessing done. Moving to training directory:
cd training/

mkdir -p models/diff_models_AIC_individual_exper_same_cell_line/
mkdir -p plots/individual_exper_same_cell_line/Cerberus plots/individual_exper_diff_cell_line/Cerberus
mkdir -p plots/individual_exper_same_cell_line/LAPA plots/individual_exper_diff_cell_line/LAPA
Rscript train_diff_models_AIC.R Cerberus
Rscript train_diff_models_AIC.R LAPA
echo "############################# 7 different logit models trained and AIC values are calculated."

# Training and predicting on the final model: logit(label~DHS + TPM + length)
for LR in Cerberus LAPA; do
	mkdir -p models/individual_exper_same_cell_line/"$LR"
	Rscript train_glm.R "$LR" # create RDS files for the model: logit(label~DHS + TPM + length)
	
	mkdir -p predicts/individual_exper_same_cell_line/"$LR"
	mkdir -p predicts/individual_exper_diff_cell_line/"$LR"
	Rscript predict_glm.R "$LR" # testing on chr2/3 LR peaks of the same cell line (K562 or GM12878) as that of the training set and creating AUC files
	Rscript predict_glm_diff_cell_line.R "$LR" # testing on chr2/3 LR peaks of the other cell line and creating AUC files
	# predictions in "predicts/" directory
	# all AUC plots in "training/plots" directory
done
echo "############################# logit(label~TPM+DHS+peak_len) trained and tested."
# SFig 6d: training/plots/Cerberus_AIC.png  training/plots/LAPA_AIC.png
# SFig 6e: training/plots/individual_exper_same_cell_line/Cerberus/GM12878_ENCSR962BVU_AUROC.pdf training/plots/individual_exper_same_cell_line/LAPA/GM12878_ENCSR962BVU_AUROC.pdf
# Making AUROC boxplot peak for cross-cell-type predictions:
Rscript plot_AUROC.R
# SFig 6f: plots/AUROC_diff_cell_line_boxplot.pdf

# Making peak statistic figures (SFig6 a, b, c)
cd .. # to main directory
mkdir -p plots
Rscript make_peak_plots.R
# SFig 6a: plots/encode_rampage_cage.pdf
# SFig 6b: plots/overlap.pdf
# SFig 6c: plots/test_train_set_peaks.pdf

echo "############################# All figure panels are created in 'plots/' and 'training/plots' directories."

