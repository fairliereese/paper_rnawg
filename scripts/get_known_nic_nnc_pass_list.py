import pandas as pd
import sys

# creates a talon pass list with only known, nic, nnc

# usage
# python known_nic_nnc_pass_list.py <abundance file> <pass list> <opref>

ab = sys.argv[1]

df = pd.read_csv(ab, sep='\t')
print(ab)
print(df.head())
df = df[['gene_ID', 'transcript_ID', 'transcript_novelty']]

pass_list = sys.argv[2]
pass_df = pd.read_csv(pass_list, header=None,
    names=['gene_ID', 'transcript_ID'])

# merge pass list in with df to get novelty of each tid
df = df.merge(pass_df, how='inner', on=['gene_ID', 'transcript_ID'])

# subset on known, nic, and nnc
novs = ['Known', 'NIC', 'NNC']
df = df.loc[df.transcript_novelty.isin(novs)]

# dump to pass list file
opref = sys.argv[3]
fname = '{}_known_nic_nnc_pass_list.csv'.format(opref)
df.to_csv(fname, header=None, index=False)
