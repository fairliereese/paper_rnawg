There is some overlap between the processing code here and in the raw data processing directory.

Notes on pre-generated files included in this directory:
* Transcription factor ids generated using [this](http://www.ensembl.org/biomart/martview/9ae34b91ac4887f7cb4e59a962bf8f87?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id_version&FILTERS=hsapiens_gene_ensembl.default.filters.go_parent_term."GO:0003700"&VISIBLEPANEL=resultspane) Biomart query.
* Human-mouse ortholog gene ids generated using [this](http://www.ensembl.org/biomart/martview/7207f9a6b715260989ef4d6aa3c1205f?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.homologs.ensembl_gene_id|hsapiens_gene_ensembl.default.homologs.ensembl_gene_id_version|hsapiens_gene_ensembl.default.homologs.ensembl_transcript_id|hsapiens_gene_ensembl.default.homologs.ensembl_transcript_id_version|hsapiens_gene_ensembl.default.homologs.mmusculus_homolog_ensembl_gene|hsapiens_gene_ensembl.default.homologs.mmusculus_homolog_associated_gene_name&FILTERS=&VISIBLEPANEL=attributepanel) Biomart query.

Before generating any of the figures, you will need to download data from the ENCODE portal and process it to generate some additional files used. Use the following Snakemake code to do so. Please note that some data processing steps in this Snakemake pipeline will download a lot of data.

```bash
snakemake \
  -s snakemake/Snakefile \
  -j 20 \
  --latency-wait 120 \
  -n
```

Currently, the figures and numbers generated for the paper are in Python notebooks roughly split by section of the paper.

Table of contents:
* [Figure 1: Data overview (S1, S2)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig1/fig1.ipynb)
* [Figure 2: Triplet feature characterization (S5, S6)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig2/fig2.ipynb)
* [Figure 3: Gene triplets (S10)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig3/fig3.ipynb)
* [Figure 4: More gene triplets and observed vs. observed major](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig4/fig4.ipynb)
* [Figure 4: MANE analyses (S11)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig4/fig_mane.ipynb)
* [Figure 5: Human-mouse comparison (S12)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig5/fig5.ipynb)
* [Figure S3: microRNA)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig_mirna/fig_mirna.ipynb)
* [Figure S7: TSS prediction](https://github.com/fairliereese/paper_rnawg/tree/master/figures/figS7)
* [Figure S8: SUPPA analyses](https://github.com/fairliereese/paper_rnawg/tree/master/figures/figS8)


Requirements:
* [Swan](https://github.com/fairliereese/paper_rnawg/tree/master/proc)
* [Cerberus](https://github.com/fairliereese/cerberus)
