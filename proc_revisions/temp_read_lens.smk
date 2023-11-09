import pandas as pd
import os
import sys
import numpy as np

def process_lr_metadata(cfg_entry, species, datasets_per_talon_run):
    """
    Concatenate metadata for each dataset from multiple species together.
    Assign each dataset to a talon run number
    """
    for i,s in enumerate(species):
        f = expand(cfg_entry, species=s)[0]
        if i == 0:
            df = pd.read_csv(f, sep='\t')
            df['species'] = s
        else:
            temp = pd.read_csv(f, sep='\t')
            temp['species'] = s
            df = pd.concat([df, temp], axis=0)

    # get number to mod each talon run by
    df['n_datasets'] = df[['dataset', 'species']].groupby('species')[['dataset']].transform('count')
    df['mod_num'] = np.ceil(df.n_datasets/datasets_per_talon_run)
    df['talon_run'] = (df.sort_values(['species', 'dataset'],
                                    ascending=[True, True])\
                                    .groupby('species')\
                                    .cumcount() + 1).to_numpy()\
                                    % df.mod_num.to_numpy()
    df['talon_run'] = df.talon_run.astype('int')
    df['max_talon_run'] = df[['talon_run', 'species']].groupby('species')[['talon_run']].transform('max')
    return df

def get_encid_from_dataset(dataset, meta, file_format):
    m = {'label_bam': 'ENCODE_alignments_id',
     'bam': 'ENCODE_unfiltered_alignments_id',
     'fastq': 'ENCODE_reads_id'}
    if file_format in list(m.keys()):
        id = meta.loc[meta.dataset == dataset, m[file_format]].values[0]
    else:
        id = meta.loc[meta.dataset == dataset, file_format].values[0]
    return id

def get_meta_df(config, species):
    meta_df = pd.DataFrame()
    for f, s in zip(list(expand(config['lr']['meta'], species=species)), species):
        temp = pd.read_csv(f, sep='\t')
        temp['species'] = s
        meta_df = pd.concat([meta_df, temp], axis=0)
    # if meta_df.dataset.duplicated.any():
    #     raise ValueError('Mouse and human dataset names not unique')
    return meta_df

configfile: 'config.yml'
datasets_per_talon_run = config['params']['talon']['datasets_per_run']
end_modes = ['tss', 'tes']
species=['mouse', 'human']
meta_df = get_meta_df(config, species)

# lr stuff
lr_df = process_lr_metadata(config['lr']['meta'],
                            species,
                            datasets_per_talon_run)

include: 'download.smk'
include: 'samtools.smk'
include: 'spike_refs.smk'
include: 'refs.smk'
include: 'talon.smk'
include: 'lapa.smk'
include: 'cerberus.smk'
include: 'gtex.smk'
include: 'encode_tss.smk'
include: 'ccre.smk'
include: 'fantom_cage.smk'
include: 'polyasite_atlas.smk'
include: 'pas.smk'
include: 'lrgasp_cage.smk'
include: 'procap.smk'


wildcard_constraints:
    dataset='|'.join([re.escape(x) for x in lr_df.dataset.tolist()]),
    species='|'.join([re.escape(x) for x in lr_df.species.unique().tolist()]),
    talon_run='|'.join([re.escape(x) for x in lr_df.talon_run.astype(str).unique().tolist()])

datasets = meta_df.loc[meta_df.species=='human'].dataset.tolist()

rule all:
    input:
        expand(config['lr']['fastq_gz'],
               species='human',
               dataset=datasets)

rule dl_lr_fastq:
    resources:
        mem_gb = 32,
        threads = 8
    params:
        encid = lambda w:get_encid_from_dataset(w.dataset,
                                                meta_df,
                                                'fastq')
    output:
        fastq = config['lr']['fastq_gz']
    shell:
        "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.fastq.gz -O {output.fastq}"
