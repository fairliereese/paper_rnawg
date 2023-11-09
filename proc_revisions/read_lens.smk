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

# TODO
datasets = meta_df.loc[meta_df.species=='human'].dataset.tolist()

################################################################################
############################## Diane's stuff ###################################
################################################################################
rule dl_lr_bam:
    resources:
        mem_gb = 32,
        threads = 8
    params:
        encid = lambda w:get_encid_from_dataset(w.dataset,
                                                meta_df,
                                                'bam')
    output:
        bam = temporary(config['lr']['bam'])
    shell:
        "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.bam -O {output.bam}"

rule dl_lr_label_bam:
    resources:
        mem_gb = 32,
        threads = 8
    params:
        encid = lambda w:get_encid_from_dataset(w.dataset,
                                                meta_df,
                                                'label_bam')
    output:
        bam = temporary(config['lr']['label_bam'])
    shell:
        "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.bam -O {output.bam}"


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

def get_col_from_meta_df(wc, col):
    temp = meta_df.copy(deep=True)
    if 'species' in wc.keys():
        temp = meta_df.loc[meta_df.species == wc['species']]
    return temp[col].tolist()

rule get_lr_read_lens:
    input:
        bams = lambda wc:expand(config['lr']['bam'],
                      species='human',
                      dataset=get_col_from_meta_df(wc, col='dataset')),
        fastqs = lambda wc:expand(config['lr']['fastq_gz'],
                        species='human',
                        dataset=get_col_from_meta_df(wc, col='dataset'))
    resources:
        mem_gb = 32,
        threads = 8
    output:
        tsv = config['lr']['read_len_meta']
    run:
        get_lr_read_lens(input.bams, input.fastqs, output.tsv)

rule all_read_lens:
    input:
        expand(rules.dl_lr_fastq.output.fastq,
               species='human',
               dataset=datasets)
