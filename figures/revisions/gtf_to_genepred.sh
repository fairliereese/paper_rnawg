gtf=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.gtf
genepred=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.genePred
biggenepred=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.bigGenePred
biggenepredsorted=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.bigGenePred_sorted
biggenepredsortedfilt=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.bigGenePred_sorted_filt
bb=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.bb
gtfToGenePred -genePredExt $gtf $genepred
genePredToBigGenePred $genepred $biggenepred
wget https://genome.ucsc.edu/goldenPath/help/examples/bigGenePred.as
sort -k1,1 -k2,2n $biggenepred > $biggenepredsorted
grep -v chrM $biggenepredsorted > $biggenepredsortedfilt
bedToBigBed -type=bed12+8 -tab -as=bigGenePred.as $biggenepredsortedfilt http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes $bb

# sample-level one
gtf=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.gtf
genepred=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.genePred
biggenepred=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.bigGenePred
biggenepredsorted=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.bigGenePred_sorted
biggenepredsortedfilt=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.bigGenePred_sorted_filt
bb=/Users/fairliereese/mortazavi_lab/data/paper_rnawg/figures/data/human/cerberus.bb
gtfToGenePred -genePredExt $gtf $genepred
genePredToBigGenePred $genepred $biggenepred
wget https://genome.ucsc.edu/goldenPath/help/examples/bigGenePred.as
sort -k1,1 -k2,2n $biggenepred > $biggenepredsorted
grep -v chrM $biggenepredsorted > $biggenepredsortedfilt
bedToBigBed -type=bed12+8 -tab -as=bigGenePred.as $biggenepredsortedfilt http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes $bb
