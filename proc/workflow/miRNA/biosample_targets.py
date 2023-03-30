import pandas as pd
import pyranges as pr


df_families = pd.read_csv(snakemake.input['family'], sep='\t')

df_counts = pd.read_csv(snakemake.input['counts'])

df_targetscan = pr.read_bed(snakemake.input['targetscan'], as_df=True)

gr_exp = pd.read_csv(snakemake.input['exp_bed'], sep='\t', header=None)
gr_exp = gr_exp.rename(columns={
    0: 'Chromosome',
    1: 'Start',
    2: 'End',
    5: 'Strand',
    9: 'cell_line',
    10: 'num_dataset'
})
gr_exp.columns = gr_exp.columns.astype(str)
gr_exp = pr.PyRanges(gr_exp)

__import__("pdb").set_trace()
