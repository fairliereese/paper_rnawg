import pandas as pd
import pyranges as pr
import cerberus
import numpy as np
from snakemake.io import expand

def get_encid_from_dataset(dataset, meta, file_format):
    m = {'label_bam': 'ENCODE_alignments_id',
     'bam': 'ENCODE_unfiltered_alignments_id',
     'fastq': 'ENCODE_reads_id'}
    if file_format in list(m.keys()):
        id = meta.loc[meta.dataset == dataset, m[file_format]].values[0]
    else:
        id = meta.loc[meta.dataset == dataset, file_format].values[0]
    return id

def format_lapa_ends(bed, o):
    df = pd.read_csv(bed, sep='\t', header=None)

    # subset on specific columns
    df = df[[0,1,2,5,7]]
    df.columns = ['Chromosome', 'Start', 'End', 'Strand', 'gene_id']

    # get stable gid
    df = cerberus.add_stable_gid(df)
    df.rename({'gene_id': 'Name'}, axis=1, inplace=True)

    # add indicator that this file came from lapa
    df['ThickStart'] = 'lapa'

    # remove regions not associated with a gene
    df = df.loc[~df.Name.str.contains('intergenic')]

    # add arbitrary unique index to name
    df['ThickEnd'] = [i for i in range(len(df.index))]
    df['Name'] = df.Name+'_'+df.ThickEnd.astype(str)
    df.drop('ThickEnd', axis=1, inplace=True)
    df = pr.PyRanges(df)
    df = df.to_bed(o, chain=True)
    return df

def get_lr_encid(wc, df):
    """
    Get the ENCID of the alignments file
    given the dataset and species name
    """
    dataset = wc.dataset
    species = wc.species
    temp = df.loc[(df.dataset==dataset)&\
                  (df.species==species)]
    return temp['ENCODE_alignments_id'].values[0]

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


def format_metadata_col(df, col, new_col):
    df[new_col] = df[col].str.lower()
    df[new_col] = df[new_col].str.replace('-', '_')
    df[new_col] = df[new_col].str.replace(' ', '_')
    return df

def process_encode_metadata(fname):
    df = pd.read_csv(fname, sep='\t')
    cols = ['Experiment accession', 'Biosample term name', 'File accession', 'Output type',
       'Biological replicate(s)', 'Technical replicate(s)']
    df = df[cols]
    df = format_metadata_col(df, 'Biosample term name', 'biosamp')
    df = format_metadata_col(df, 'Output type', 'output')

    # get biorep number for each experiment
    keep_cols = ['Experiment accession', 'biosamp']
    df.sort_values(by='Experiment accession', ascending=True)
    temp = df[keep_cols].drop_duplicates(keep='first')
    temp['biorep'] = temp[keep_cols].groupby('biosamp').cumcount()+1
    temp.drop('biosamp', axis=1, inplace=True)

    # merge in biorep
    df = df.merge(temp, how='left', on='Experiment accession')

    return df

def add_feat(df, col, kind, as_index=False, as_number=False,
             drop_gid=True,
             drop_triplet=True):
    """
    Add ic, tss, tes info to a df from the transcript id

    Parameters:
        df (pandas DataFrame): Df w/ tid in there
        col (str): Which column to use
        kind (str): {'tss', 'ic', 'tes'}
        as_index (bool): Whether to replace current index with
            new feat
        as_number (bool): Just add the number of element, not
            geneid_#
        drop_gid (bool): Drop added gene id
        drop_triplet (bool): Drop added triplet

    Returns:
        df (pandas DataFrame): DF w/ additional "kind" col
    """

    if col == 'index':
        df['temp_tid'] = df.index.tolist()
    else:
        df['temp_tid'] = df[col].tolist()
    df['triplet'] = df.temp_tid.str.split('[', expand=True)[1].str.split(']', expand=True)[0]
    df['temp_gid'] = df.temp_tid.str.split('[', expand=True)[0]
    if kind == 'tss':
        ind = 0
    elif kind == 'ic':
        ind = 1
    elif kind == 'tes':
        ind = 2
    df[kind] = df.triplet.str.split(',', expand=True)[ind]
    if not as_number:
        df[kind] = df.temp_gid+'_'+df[kind]
    drop_cols = list(['temp_tid'])
    if drop_gid:
        drop_cols += ['temp_gid']
    if drop_triplet:
        drop_cols += ['triplet']
    df.drop(drop_cols, axis=1, inplace=True)

    if as_index:
        df.reset_index(inplace=True, drop=True)
        df.index = df[kind]
        df.drop(kind, axis=1, inplace=True)
    return df

def get_det_table(meta_file,
                  enc_meta_file,
                  filt_ab,
                  min_tpm,
                  ofile):

    # get dataset<->experiment<->encode biosample shorthand
    meta = pd.read_csv(meta_file, sep='\t')
    meta = meta[['ENCODE_experiment_id', 'dataset']].rename({'ENCODE_experiment_id': 'Experiment accession'}, axis=1)

    enc_meta = process_encode_metadata(enc_meta_file)
    enc_meta = enc_meta[['Experiment accession', 'biosamp', 'biorep']].drop_duplicates()

    meta = meta.merge(enc_meta, how='left', on='Experiment accession')
    meta = meta[['dataset', 'biosamp', 'biorep']]
    meta.set_index('dataset', inplace=True)

    # get counts / tss
    df = pd.read_csv(filt_ab, sep='\t')
    df = add_feat(df,
                  'annot_transcript_id',
                  kind='tss',
                  as_index=True)
    drop_cols = ['gene_ID', 'transcript_ID',
         'annot_gene_id', 'annot_transcript_id',
         'annot_gene_name', 'annot_transcript_name',
         'n_exons', 'length', 'gene_novelty',
         'transcript_novelty', 'ISM_subtype']
    df.drop(drop_cols, axis=1, inplace=True)
    df = df.reset_index()
    df = df.groupby('tss').sum()

    # get the tpm
    df = df.transpose()
    df['total'] = df.sum(axis=1)
    df = df.divide(df.total, axis=0)*1e6
    df.drop('total', axis=1, inplace=True)

    # enforce min tpm
    df = df>=min_tpm

    # add metadata so we can merge on biosample as dictated by encode
    df = df.merge(meta, how='left', left_index=True, right_index=True)
    df.reset_index(drop=True, inplace=True)
    df = df.groupby(['biosamp', 'biorep']).max()

    df.to_csv(ofile, sep='\t')

def get_lr_tss(ca_h5,
               det_mat,
               biosamp,
               biorep,
               ofile):

    # get tsss
    ca = cerberus.read(ca_h5)
    tss = ca.tss.copy(deep=True)

    # get det info
    df = pd.read_csv(det_mat, sep='\t')
    df = df.loc[(df.biosamp==biosamp)&(df.biorep.astype(int)==int(biorep))]
    df = df.drop(['biosamp', 'biorep'], axis=1).transpose()
    df.columns = ['temp']
    df = df.loc[df['temp']==True]
    tss_ids = df.index.tolist()

    # get bed file of det tsss
    tss = tss.loc[tss.Name.isin(tss_ids)]
    tss = pr.PyRanges(tss)
    tss.to_bed(ofile)

def intersect_ccre(ccre_file,
                   bed_files,
                   ofile):
    ccre = pr.read_bed(ccre_file)

    # get concatenation of all bed files
    df = pd.DataFrame()
    for f in bed_files:
        temp = pd.read_csv(f, sep='\t')
        df = pd.concat([df, temp], axis=0)

    # merge w/ ccre data
    df = pr.PyRanges(df)

    ccre = ccre.join(df,
         how='left',
         slack=0,
         suffix='_other').df

    # get support for each ccre by each assay
    temp = ccre[['Name', 'assay']].drop_duplicates()
    temp = temp.sort_values(by='assay', ascending=True)
    temp = temp.loc[temp.assay!='-1']
    assay_list = [temp.assay.unique()]+[np.nan]
    temp = temp.groupby('Name').agg({'assay':','.join}, axis=0).reset_index()
    temp['assay'] = temp['assay']+',ccre'

    # now make table w/ support for each ccre
    ccre = pr.read_bed(ccre_file).df
    ccre = ccre.merge(temp, how='left', on='Name')
    ccre.loc[ccre.assay.isnull(), 'assay'] = 'ccre'

    # and save
    ccre.to_csv(ofile, sep='\t', index=False)

def filter_gtex_gtf(gtf, oname):
    df = pr.read_gtf(gtf, as_df=True)
    df = df.loc[~df.gene_id.str.contains('chr')]
    # df = df.loc[df.gene_id.str.contains('ENSG')]
    # inds = df.loc[df.gene_id.str.contains('chr')].index.tolist()
    # df.loc[inds, 'gene_id'] = df.loc[inds, 'gene_id'].str.split('_', n=1, expand=True)[1]
    df = pr.PyRanges(df)
    df.to_gtf(oname)
