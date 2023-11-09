import pandas as pd
import os
import sys
import numpy as np



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
species=['mouse', 'human']
meta_df = get_meta_df(config, species)

datasets = meta_df.loc[meta_df.species=='human'].dataset.tolist()

rule all:
    input:
        expand(rules.dl_lr_fastq.output.fastq,
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
