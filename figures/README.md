Before generating any of the figures, you will need to download data from the ENCODE portal and process it to generate some additional files used.

There is some overlap between the processing code here and in the raw data processing directory.

Notes on pre-generated files included in this directory:
* Biomart TF ids generated using [this](http://www.ensembl.org/biomart/martview/9ae34b91ac4887f7cb4e59a962bf8f87?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id_version&FILTERS=hsapiens_gene_ensembl.default.filters.go_parent_term."GO:0003700"&VISIBLEPANEL=resultspane) query.

```bash
snakemake \
  -s snakemake/human_snakefile.smk \
  -j 10 \
  --latency-wait 120 \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END, --time=72:00:00" -n

snakemake \
  -s snakemake/Snakefile \
  -j 10 \
  --latency-wait 120 \
  -n
```
