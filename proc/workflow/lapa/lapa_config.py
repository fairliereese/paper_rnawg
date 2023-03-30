import pandas as pd


df = snakemake.params.df_lr[['sample']].rename(
    columns={'sample': 'dataset'})
df = df.reset_index().rename(
    columns={'File accession': 'sample'})

bam_path = snakemake.config['lr']['bam']
df['path'] = df['sample'].map(lambda x: bam_path.replace('{encode_id}', x))

df.to_csv(snakemake.output['config'], index=False)
