## Machine learning models predict the support for long-read TSS peaks by other TSS-annotating assays and in a cross-cell type manner

This section contains the code to train/test logistic regression models to predict the overlap (support) for a peak in different TSS assays.\
The steps to perform preprocessing, training and testing the logit models are in `run.sh`. `.sh` files and `bigWigAverageOverBed` must have execution permission. \
`./run.sh` should generate all the panels in SFig 6. \
You will need some libraries and tools including  PRROC, dplyr, ggplot2, RColorBrewer in R and bedtools on your machine.

After running this script, Supp. Fig. 6 panels should be created in the following directories: \
SFig 6a: plots/encode_rampage_cage.pdf \
SFig 6b: plots/overlap.pdf \
SFig 6c: plots/test_train_set_peaks.pdf \
SFig 6d: training/plots/Cerberus_AIC.png  training/plots/LAPA_AIC.png \
SFig 6e: training/plots/individual_exper_same_cell_line/Cerberus/GM12878_ENCSR962BVU_AUROC.pdf training/plots/individual_exper_same_cell_line/LAPA/GM12878_ENCSR962BVU_AUROC.pdf \
SFig 6f: plots/AUROC_diff_cell_line_boxplot.pdf
