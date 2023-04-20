## Input files
1. `cerberus_transcript_novelty.tsv`
2. `gtex_to_cerberus_id.tsv`
3. `suppa/localEvents/*.ioe` and `suppa/localEvents_GTEx/*.ioe` from the SUPPA output (see S7A).

## Sankey plot
Run `Rscript sankey_plot.R` to generate `sankey_plot.pdf` showing the number of transcripts classified as Known, Novel In Catalog (NIC), Novel Not in Catalog (NNC), Unspliced, or Missing by Cerberus and GTEx.
