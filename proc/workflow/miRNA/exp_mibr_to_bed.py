import pdb
import pandas as pd
import pyranges as pr


df = pd.read_csv(snakemake.input['tsv'], sep='\t') \
       .rename(columns={
           'CHROM': 'Chromosome',
           'start_Ref38': 'Start',
           'end_Ref38': 'End',
           'STRAND': 'Strand'
       })

gr = pr.PyRanges(df[[
    'Chromosome', 'Start', 'End', 'MIRNA', 'NDATASET', 'Strand', 'DATASET',
    'NLINES', 'LINES', 'REGION', 'GENES'
]])
gr.to_bed(snakemake.output['all_bed'])

gr[~gr.MIRNA.isna()].to_bed(snakemake.output['family_bed'])
