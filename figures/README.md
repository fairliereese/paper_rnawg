There is some overlap between the processing code here and in the raw data processing directory.

Notes on pre-generated files included in this directory:
* Transcription factor ids generated using [this](http://www.ensembl.org/biomart/martview/9ae34b91ac4887f7cb4e59a962bf8f87?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id_version|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id_version&FILTERS=hsapiens_gene_ensembl.default.filters.go_parent_term."GO:0003700"&VISIBLEPANEL=resultspane) Biomart query.
* Human-mouse ortholog gene ids generated using [this](http://www.ensembl.org/biomart/martview/7207f9a6b715260989ef4d6aa3c1205f?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.homologs.ensembl_gene_id|hsapiens_gene_ensembl.default.homologs.ensembl_gene_id_version|hsapiens_gene_ensembl.default.homologs.ensembl_transcript_id|hsapiens_gene_ensembl.default.homologs.ensembl_transcript_id_version|hsapiens_gene_ensembl.default.homologs.mmusculus_homolog_ensembl_gene|hsapiens_gene_ensembl.default.homologs.mmusculus_homolog_associated_gene_name&FILTERS=&VISIBLEPANEL=attributepanel) Biomart query.

Before generating any of the figures, you will need to download data from the ENCODE portal and process it to generate some additional files used. Use the following Snakemake code to do so. Please note that some data processing steps in this Snakemake pipeline will download a lot of data.
<!--

```python
import os
import sys
from encoded_client.encoded import ENCODED

p = os.path.dirname(os.getcwd())
sys.path.append(p)

from scripts.utils import *   
df = get_lr_exp_meta('human')
df['species'] = 'human'
df2 = get_lr_exp_meta('mouse')
df2['species'] = 'mouse'
df = pd.concat([df, df2], axis=0)
output_types = ['unfiltered alignments', 'alignments', 'reads']
df = df.loc[df.output_type.isin(output_types)]
df.to_csv('ref/lr_file_ids.tsv', sep='\t', index=False)
``` -->

```bash
snakemake \
  -s snakemake/Snakefile \
  -j 30 \
  --latency-wait 120 \
  --use-conda \
  -n

conda activate viz_snakemake
snakemake -s snakemake/Snakefile --dag | dot -Tpng > ruledag.png

snakemake \
  -s snakemake/Snakefile \
  -j 60 \
  --latency-wait 120 \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n
```

Currently, the figures and numbers generated for the paper are in Python notebooks roughly split by section of the paper.

Table of contents:
* [Figure 1: Data overview (S1, S2)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig1/fig1.ipynb)
* [Figure 2: Triplet feature characterization (S5, S6)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig2/fig2.ipynb)
* [Figure 3: Gene triplets (S10)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig3/fig3.ipynb)
* [Figure 4: More gene triplets and observed vs. observed major](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig4/fig4.ipynb)
* [Figure 4: MANE analyses (S11)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig4/fig_mane.ipynb)
* [Figure 5: Human-mouse comparison (S12)](https://github.com/fairliereese/paper_rnawg/blob/master/figures/fig5/fig5.ipynb)
* [Figure S2: microRNA](https://github.com/fairliereese/paper_rnawg/blob/master/figures/figS2/fig_mirna.ipynb)
* [Figure S6: TSS prediction](https://github.com/fairliereese/paper_rnawg/tree/master/figures/figS6)
* [Figure S8: SUPPA analyses](https://github.com/fairliereese/paper_rnawg/tree/master/figures/figS8/figS8.ipynb)
* [Protein prediction](https://github.com/fairliereese/paper_rnawg/blob/master/figures/protein_pred/pp_overview.ipynb)


Requirements:
* [Swan](https://github.com/fairliereese/swan_vis)
* [Cerberus](https://github.com/fairliereese/cerberus)
* [ENCODED client](https://github.com/detrout/encoded_client)
