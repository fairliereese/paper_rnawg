```bash
as_file=track_meta.as
ifile=test_meta.bgp
chrom_sizes=ref/human/chrom_sizes.tsv
ofile=data/human/tracks/huvec_test_new_as.bb
bedToBigBed -type=bed12+8 -tab -as=$as_file $ifile $chrom_sizes $ofile
```
