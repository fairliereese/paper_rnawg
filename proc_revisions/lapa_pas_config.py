import pandas as pd


df = snakemake.params.df_pas[['Biosample term name']].rename(
    columns={'Biosample term name': 'dataset'})
df = df.reset_index().rename(
    columns={'File accession': 'sample'})
df = df.loc[df.species==snakemake.wildcards.species]

bam_path = snakemake.config['pas']['bam']
df['path'] = df['sample'].map(lambda x: bam_path.replace('{encid}', x))

df.to_csv(snakemake.output['config'], index=False)