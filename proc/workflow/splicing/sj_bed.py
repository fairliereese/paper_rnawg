import pyranges as pr
import pandas as pd


df = pd.read_csv(snakemake.input['ssj'], sep='\t', header=None)
intervals = df[0].str.split('_')

df = pd.DataFrame({
    'Chromosome': intervals.map(lambda x: ''.join(x[:-3])),
    'Start': intervals.map(lambda x: x[-3]),
    'End': intervals.map(lambda x: x[-2]),
    'Name': '.',
    'Score': df[3],
    'Strand': intervals.map(lambda x: x[-1])
})
df['Start'] = df['Start'].astype('int') - 1

pr.PyRanges(df).to_bed(snakemake.output['sj'])
