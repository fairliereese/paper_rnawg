mkdir data
cd data

# download
xargs -L 1 curl -O -J -L < ../files.txt

mv metadata.tsv ../
cd ../
mkdir plots

# process and make figures
Rscript miRNA_figure.R
