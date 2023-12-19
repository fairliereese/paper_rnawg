import pandas as pd
import numpy as np
import scipy.stats as st
import seaborn as sns
import sys
import os
import gseapy as gp
import matplotlib.pyplot as plt
import swan_vis as swan
import yaml
from snakemake.io import expand


def read_annot_to_counts(df):
    """
    Convert a long-form (one row = one read)
    read_annot style dataframe to a counts matrix 
    """
    # compute transcript counts from the sampled df
    temp = df.copy(deep=True)
    temp = temp.groupby(['dataset', 'annot_transcript_id']).count().reset_index()
    temp.set_index('annot_transcript_id', inplace=True)
    temp = temp.pivot(columns='dataset', values='read_name') 
    temp.columns.name = ''
    temp.fillna(0, inplace=True)
    temp.reset_index(inplace=True)
    return temp

def add_ab_metdata(df, 
                   ab):
    """
    Add in gene and transcript metadata from a talon-style 
    abundance matrix
    """
    # merge in the transcript / gene info from the abundance file
    ab_df = pd.read_csv(ab, sep='\t')
    non_dataset_columns = ['gene_ID', 'transcript_ID', 'annot_gene_id',
                       'annot_transcript_id', 'annot_gene_name',
                       'annot_transcript_name', 'n_exons', 'length',
                       'gene_novelty', 'transcript_novelty', 'ISM_subtype']
    ab_df = ab_df[non_dataset_columns]
    df = df.merge(ab_df, how='left', on='annot_transcript_id')
    return df
    

def subsample_ab(talon_ab, 
                       depth, 
                       ofile):
    """
    Subsample WTC11 data to a specified depth from the abundance matrix
    and output a TALON abundance file.
    """
    
    df = pd.read_csv(talon_ab, sep='\t')    
    
    # limit to just wtc11  datasets
    wtc11_cols = ['wtc11_1_3', 'wtc11_1_2', 'wtc11_1_1']
    df = df[['annot_transcript_id']+wtc11_cols]
    
    counts1 = df.set_index('annot_transcript_id').sum(axis=1).sum(axis=0) 
    
    # melt to get counts for each transcript, dataset pair
    df = pd.melt(df, id_vars='annot_transcript_id',
                 var_name='dataset',
                 value_name='counts')
    
    # get a df w/ one entry for transcript, dataset for each count
    reads = []
    tids = []
    datasets = []
    for ind, entry in df.iterrows():
        counts = entry.counts

        tid = [entry.annot_transcript_id for i in range(counts)]
        dataset = [entry.dataset for i in range(counts)]

        tids += tid
        datasets += dataset

    df = pd.DataFrame()
    df['annot_transcript_id'] = tids
    df['dataset'] = datasets
    counts2 = len(df.index)
    assert counts1 == counts2
    
    # sample to target depth
    sample_df = df.sample(frac=depth,
                          replace=False)
    sample_df['read_name'] = sample_df.index  
    
    temp = read_annot_to_counts(sample_df)
    temp = add_ab_metdata(temp, talon_ab)
    
    temp.to_csv(ofile, sep='\t', index=False)
    
def subsample_read_annot(read_annot,
                         talon_ab,
                         depth,
                         ofile):
    """
    Subsample WTC11 data to a specified depth from the read_annot file
    and output a TALON abundance file.
    """
    
    # read in the transcript id, dataset, and the read name from the read annot
    df = pd.read_csv(read_annot, sep='\t', 
                     usecols=[0,1,12])
    
    # limit to just the wtc11 datasets
    wtc11_cols = ['wtc11_1_3', 'wtc11_1_2', 'wtc11_1_1']
    df = df.loc[df.dataset.isin(wtc11_cols)]
    
    # sample to target depth
    sample_df = df.sample(frac=depth,
                          replace=False)
    
    temp = read_annot_to_counts(sample_df)
    temp = add_ab_metdata(temp, talon_ab)
    
    temp.to_csv(ofile, sep='\t', index=False)
    