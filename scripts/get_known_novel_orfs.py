from pyfaidx import Fasta
import pandas as pd
import sys

# usage python get_known_novel_orfs.py tama_orfs.fa gencode_pc_translations.fa opref

pred_orfs = sys.argv[1]
known_orfs = sys.argv[2]
opref = sys.argv[3]
fname = opref+'_known_novel_orf.fa'

def make_fa_df(orfs):
    beeps = Fasta(orfs)  
    
    keys = list(beeps.keys())
    items = list(beeps.values())
    items = [str(item) for item in items]
    
    df = pd.DataFrame()
    df['id'] = keys
    df['seq'] = items
    
    return df

# novel orfs
df = make_fa_df(pred_orfs)

df['novel'] = df.id.str.contains('ENCODEHT')
df[['id', 'novel']].groupby('novel').count()

# drop known bois
df = df.loc[df.novel == True]
df['id'] = df.id.str.split(';', expand=True)[1]
df['id'] = '>'+df.id
df.drop('novel', axis=1, inplace=True)

ids = df.id.to_frame()
ids.columns = ['entry']
ids['type'] = 0
seqs = df.seq.to_frame()
seqs.columns = ['entry']
seqs['type'] = 1

df = pd.concat([ids, seqs], axis=0)

df['sort_col'] = df.index.tolist()
df.sort_values(by=['sort_col', 'type'], inplace=True, ascending=[True, True])
df.drop(['sort_col', 'type'], axis=1, inplace=True)

df.to_csv(fname, index=False, header=None)

# known orfs
df = make_fa_df(known_orfs)

df['id'] = df.id.str.split('|', expand=True)[1]
df['id'] = '>'+df.id
ids = df.id.to_frame()
ids.columns = ['entry']
ids['type'] = 0
seqs = df.seq.to_frame()
seqs.columns = ['entry']
seqs['type'] = 1

df = pd.concat([ids, seqs], axis=0)

df['sort_col'] = df.index.tolist()
df.sort_values(by=['sort_col', 'type'], inplace=True, ascending=[True, True])
df.drop(['sort_col', 'type'], axis=1, inplace=True)

df.to_csv(fname, index=False, header=None, mode='a')