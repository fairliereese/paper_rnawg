import pandas as pd
import pyranges as pr


df = pd.concat([
    pr.read_bed(i, as_df=True)
    for i in snakemake.input['sj']
])

df = df.groupby(['Chromosome', 'Start', 'End', 'Strand'], observed=True) \
    .sum().reset_index()
df['Name'] = '.'

pr.PyRanges(df).to_bed(snakemake.output['sj'])
