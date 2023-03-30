import pandas as pd

df_mirna = snakemake.params['df_mirna']

df_count = pd.read_csv(snakemake.input['counts'], sep='\t')

df = pd.DataFrame()
df['gene_id'] = df_count['GENCODE_ID']
df['gene_name'] = df_count['GENCODE_Name']

for file_acc, row in df_mirna.iterrows():
    biosample = row['Biosample.term.name']

    if row['Biosample.term.name'] in df.columns:
        df[biosample] += df_count[file_acc]
    else:
        df[biosample] = df_count[file_acc]

df.to_csv(snakemake.output['counts'], index=False)
