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
import pyranges as pr

p = os.getcwd()
sys.path.append(p)

# from .utils import *
from utils import *

def add_ss_type_to_intron(df):
    """
    Given a bed-style df, add "ss_3" and "ss_5" columns to
    indicate which splice site each coordinate is used

    Parameters:
        df (pandas DataFrame): Bed-style DF of introns

    Returns:
        df (pandas DataFrame): Bed-style DF of introns w/ ss_3 and ss_5 columns added
    """

    df['ss_5'] = np.nan
    df.loc[df.Strand=='+', 'ss_5'] = df.loc[df.Strand=='+', ['Start', 'End']].min(axis=1)
    df.loc[df.Strand=='-', 'ss_5'] = df.loc[df.Strand=='-', ['Start', 'End']].max(axis=1)

    df['ss_3'] = np.nan
    df.loc[df.Strand=='+', 'ss_3'] = df.loc[df.Strand=='+', ['Start', 'End']].max(axis=1)
    df.loc[df.Strand=='-', 'ss_3'] = df.loc[df.Strand=='-', ['Start', 'End']].min(axis=1)

    assert len(df.loc[(df.ss_3<df.ss_5)&(df.Strand=='+')].index) == 0
    assert len(df.loc[(df.ss_3>df.ss_5)&(df.Strand=='-')].index) == 0

    return df


def intron_to_ss(df, id_cols=None):
    """
    Get splice site coordinates from intron coordinates in bed format

    Parameters:
        df (pandas DataFrame): Pandas DF of intron coordinates in bed format
        id_cols (None or list of str): List of columns to use as ss identifier
            during melt, otherwise None

    Returns:
        df (pandas DataFrame): Pandas DF of splice site coordinates in semi
            bed format (ie no end positions, just starts, as these are single bp)
    """


    # since these are intron coords, the start defines a 3' ss
    # and the end defines a 5' ss.
    # df.rename({'Start':'ss_3', 'End':'ss_5'}, axis=1, inplace=True)
#     df['ss_5'] = np.nan
#     df.loc[df.Strand=='+', 'ss_5'] = df.loc[df.Strand=='+', ['Start', 'End']].min(axis=1)
#     df.loc[df.Strand=='-', 'ss_5'] = df.loc[df.Strand=='-', ['Start', 'End']].max(axis=1)

#     df['ss_3'] = np.nan
#     df.loc[df.Strand=='+', 'ss_3'] = df.loc[df.Strand=='+', ['Start', 'End']].max(axis=1)
#     df.loc[df.Strand=='-', 'ss_3'] = df.loc[df.Strand=='-', ['Start', 'End']].min(axis=1)

    df = add_ss_type_to_intron(df)

    assert len(df.loc[(df.ss_3<df.ss_5)&(df.Strand=='+')].index) == 0
    assert len(df.loc[(df.ss_3>df.ss_5)&(df.Strand=='-')].index) == 0

    df.drop(['Start', 'End'], axis=1, inplace=True)

    if id_cols:
        id_cols += ['Chromosome', 'Strand']
    else:
        id_cols = ['Chromosome', 'Strand']

    df = df.melt(id_vars=id_cols,
                 var_name='ss_type',
                 value_name='Start')

    # remove duplicates, which would result from the same
    # ss being used in different sjs
    df = df.drop_duplicates()

    return df

def get_source_table(df):
    """
    Get a melted form table for each entry in a tss, ic, or tes table
    for each form of support for each triplet feature.

    Parameters:
        df (pandas DataFrame): DataFrame of tsss, ics, or tess

    Returns:
        df (pandas DataFrame): Long-form DataFrame of support for each tss, ic, or tes
    """
    keep_cols = ['Name', 'source']
    df = df[keep_cols].copy(deep=True)
    df['list_source'] = df.source.str.split(',')
    df = df.explode('list_source')
    df.drop('source', axis=1, inplace=True)

    return df

# chatgpt wrote this for me thanx chatgpt
def sequential_pairs(x):
    """
    Get sequential pairs of tuples in list.
    Example: [1,2,3,4] -> [(1,2),(3,4)]
    """
    p = []
    for i in range(0, len(x) - 1, 2):
        p.append((x[i], x[i + 1]))
    return p

def explode_ic(ic):
    """
    Explode an ic df to long form with splice junction entries
    """
    # remove the monoexonic entries
    ic = ic.loc[~(ic.Coordinates == '-')]

    # explode into series of ss coords
    keep_cols = ['Chromosome', 'Coordinates',
                 'Strand', 'gene_id',
                 'Name']
    df = ic.copy(deep=True)
    df = df[keep_cols]
    df['ss_coords'] = df.Coordinates.str.split('-')

    # get pairs of sss to form sjs
    df['sj_coords'] = df.ss_coords.apply(sequential_pairs)
    df = df.explode('sj_coords')
    df.drop(['Coordinates', 'ss_coords'], axis=1, inplace=True)

    return df

def get_ss_sj_from_ic(ic, ref_sources, how):
    ic = ic.copy(deep=True)

    # get coords of each splice site in each splice junction
    df = explode_ic(ic)
    df['Start'] = df['sj_coords'].str[0].astype(int)
    df['End'] = df['sj_coords'].str[1].astype(int)
    df.drop('sj_coords', axis=1, inplace=True)

    # label sss as 5' or 3' and melt
    if how == 'ss':
        # assert len(df.loc[(df.Start>df.End)&(df.Strand=='+')].index) == 0
        # # since these are intron coords, the start defines a 3' ss
        # # and the end defines a 5' ss
        # df.rename({'Start':'ss_3', 'End':'ss_5'}, axis=1, inplace=True)
        # id_cols = ['Chromosome', 'Strand', 'gene_id', 'Name']
        # df = df.melt(id_vars=id_cols,
        #              var_name='ss_type',
        #              value_name='Start')
        df = intron_to_ss(df, ['gene_id', 'Name'])

    # for sjs, reorder according to min and max coords
    # in bed standard format
    elif how == 'sj':
        df['temp_Start'] = df.Start
        df['temp_End'] = df.End
        df['Start'] = df[['temp_Start', 'temp_End']].min(axis=1)
        df['End'] = df[['temp_Start', 'temp_End']].max(axis=1)
        df.drop(['temp_Start', 'temp_End'], axis=1, inplace=True)

    # df to hold ic to ss or sj info
    ic_df = df.copy(deep=True)

    # merge source info in w/ coord info
    df2 = get_source_table(ic)
    df = df.merge(df2, how='left', on=['Name'])

    # figure out novelty and source of each ss / sj
    df.drop('Name', axis=1, inplace=True)
    df.drop_duplicates(inplace=True)
    gb_cols = ['Chromosome', 'Strand', 'gene_id', 'Start']
    if how == 'ss':
        gb_cols += ['ss_type']
    elif how == 'sj':
        gb_cols += ['End']
    df.rename({'list_source': 'source'},
              axis=1, inplace=True)
    df['novelty'] = df.source.isin(ref_sources).map({True: 'Known',
                                                     False: 'Novel'})
    df = df.groupby(gb_cols).agg(','.join).reset_index()
    df = cerberus.update_novelty(df)

    # add novelty and support information to the ic / (ss or sj) df
    merge_cols = ['Chromosome', 'Strand', 'gene_id', 'Start']
    if how == 'ss':
        merge_cols += ['ss_type']
    elif how == 'sj':
        merge_cols += ['End']
    ic_df = ic_df.merge(df, how='left',
                        on=merge_cols)


    return df, ic_df

def get_sj_from_ic(ic, ref_sources):
    """
    Get a splice junction table from an intron chain table.
    Retain source and novelty information.

    Parameters:
        ic (pandas DataFrame): DataFrame formatted as cerberus ic table
        ref_sources (list of str): List of sources to use as references

    Returns:
        df (pandas DataFrame): DataFrame with entries for each splice junction
        ic_df (pandas DataFrame): DataFrame with entries for each splice junction /
            intron chain combination

    """
    return get_ss_sj_from_ic(ic, ref_sources, 'sj')

def get_ss_from_ic(ic, ref_sources):
    """
    Get a splice site table from an intron chain table.
    Retain source and novelty information.

    Parameters:
        ic (pandas DataFrame): DataFrame formatted as cerberus ic table
        ref_sources (list of str): List of sources to use as references

    Returns:
        df (pandas DataFrame): DataFrame with entries for each splice site
        ic_df (pandas DataFrame): DataFrame with entries for each splice site /
            intron chain combination
    """
    return get_ss_sj_from_ic(ic, ref_sources, 'ss')

def get_fusion_sj_table(ab, ics, ref_ics, ver, include_novel=True):
    gtf_df, _, _ = get_gtf_info(ver=ver, how='gene', add_stable_gid=True)
    gtf_df = gtf_df[['gid_stable', 'gname']]
    gtf_df.head()

    talon_df = pd.read_csv(talon_filt_ab, sep='\t')
    talon_df['gid'] = cerberus.get_stable_gid(talon_df, 'annot_gene_id')

    fusion_gids = talon_df.loc[talon_df.gene_novelty=='Fusion', 'gid'].tolist()
    known_gids = talon_df.loc[talon_df.gene_novelty=='Known', 'gid'].tolist()

    df = pd.read_csv(ics, sep='\t')
    df['gene_id'] = df.Name.str.split('_', expand=True)[0]
    df['source'] = 'lapa'

    fusion_df = df.loc[df.gene_id.isin(fusion_gids)].copy(deep=True)
    # known_df = df.loc[df.gene_id.isin(known_gids)]
    known_df = pd.read_csv(ref_ics, sep='\t')
    known_df['gene_id'] = known_df.Name.str.split('_', expand=True)[0]
    known_df['source'] = 'lapa'

    _, f_sj_ic_df = get_sj_from_ic(fusion_df, ['lapa'])
    f_sj_ic_df.drop(['source', 'novelty'], axis=1, inplace=True)

    _, k_sj_ic_df = get_sj_from_ic(known_df, ['lapa'])
    k_sj_ic_df.drop(['source', 'novelty'], axis=1, inplace=True)

    f_sj_ic_df = f_sj_ic_df.merge(k_sj_ic_df, how='left', on=['Chromosome', 'Strand', 'Start', 'End'],
                              suffixes=('', '_known'))
    f_sj_ic_df['known'] = False
    f_sj_ic_df.loc[f_sj_ic_df.Name_known.notnull(), 'known'] = True

    # total # splice sjs / ic
    temp = f_sj_ic_df[['Name', 'Chromosome', 'Strand', 'Start', 'End']].copy(deep=True)
    temp['sj'] = temp.Name+temp.Chromosome+temp.Strand+temp.Start.astype(str)+temp.End.astype(str)
    temp = temp[['Name', 'sj']]
    temp = temp.groupby('Name').nunique().reset_index()
    temp.rename({'sj':'n_total_sj'}, axis=1, inplace=True)
    temp.head()

    # splice junctions / ic from each gene
    temp2 = f_sj_ic_df[['Name', 'known', 'Start', 'End', 'gene_id_known']].copy(deep=True)
    temp2.drop_duplicates(inplace=True)
    temp2.drop('End', axis=1, inplace=True)
    if include_novel == True:
        dropna = False
    else:
        dropna = True
    temp2 = temp2.groupby(['Name', 'known', 'gene_id_known'], dropna=dropna).count().reset_index()
    temp2.rename({'gene_id_known': 'gene_id', 'Start': 'n_sj'}, axis=1, inplace=True)
    temp2.head()

    # number each gene as 1' 2' or 3'
    temp2 = temp2.sort_values(by='n_sj', ascending=False)
    temp2['rank'] = temp2.sort_values(by='n_sj', ascending=False).groupby(['Name']).cumcount()+1
    temp2.head()

    # number of genes w/ intersecting sjs / ic
    temp3 = temp2[['Name', 'gene_id']].groupby('Name', dropna=dropna).nunique(dropna=dropna).reset_index()
    temp3.rename({'gene_id': 'n_genes'}, axis=1, inplace=True)
    temp3.head()

    temp5 = temp.merge(temp2, on='Name')
    temp5 = temp5.merge(temp3, on='Name')
    temp5['perc'] = (temp5['n_sj']/temp5['n_total_sj'])*100
    return f_sj_ic_df, temp5

def get_fusion_ss_table(ab, gtf, ref_ics, ver, include_novel=True):
    gtf_df, _, _ = get_gtf_info(ver=ver, how='gene', add_stable_gid=True)
    gtf_df = gtf_df[['gid_stable', 'gname']]
    gtf_df.head()

    talon_df = pd.read_csv(ab, sep='\t')
    talon_df['gid'] = cerberus.get_stable_gid(talon_df, 'annot_gene_id')

    fusion_gids = talon_df.loc[talon_df.gene_novelty=='Fusion', 'gid'].tolist()
    known_gids = talon_df.loc[talon_df.gene_novelty=='Known', 'gid'].tolist()

    # df = pd.read_csv(ics, sep='\t')
    # df['gene_id'] = df.Name.str.split('_', expand=True)[0]
    # df['source'] = 'lapa'
    df = pr.read_gtf(gtf, rename_attr=True, duplicate_attr=True)
    df = cerberus.get_ic(df)
    df.rename({'transcript_id':'Name', 'ic': 'Coordinates'}, axis=1, inplace=True)
    df['source'] = 'lapa'


    fusion_df = df.loc[df.gene_id.isin(fusion_gids)].copy(deep=True)
    known_df = pd.read_csv(ref_ics, sep='\t')
    known_df['gene_id'] = known_df.Name.str.split('_', expand=True)[0]
    known_df['source'] = 'lapa'

    _, f_ss_ic_df = get_ss_from_ic(fusion_df, ['lapa'])
    f_ss_ic_df.drop(['source', 'novelty'], axis=1, inplace=True)

    _, k_ss_ic_df = get_ss_from_ic(known_df, ['lapa'])
    k_ss_ic_df.drop(['source', 'novelty'], axis=1, inplace=True)

    f_ss_ic_df = f_ss_ic_df.merge(k_ss_ic_df, how='left', on=['Chromosome', 'Strand', 'Start', 'ss_type'],
                              suffixes=('', '_known'))
    f_ss_ic_df['known'] = False
    f_ss_ic_df.loc[f_ss_ic_df.Name_known.notnull(), 'known'] = True

    # total # splice sss / ic
    temp = f_ss_ic_df[['Name', 'Chromosome', 'Strand', 'Start']].copy(deep=True)
    temp['ss'] = temp.Name+temp.Chromosome+temp.Strand+temp.Start.astype(str)
    temp = temp[['Name', 'ss']]
    temp = temp.groupby('Name').nunique().reset_index()
    temp.rename({'ss':'n_total_ss'}, axis=1, inplace=True)
    temp.head()

    # splice sites / ic from each gene
    temp2 = f_ss_ic_df[['Name', 'known', 'Start', 'gene_id_known']].copy(deep=True)
    temp2.drop_duplicates(inplace=True)
    if include_novel == True:
        dropna = False
    else:
        dropna = True
    temp2 = temp2.groupby(['Name', 'known', 'gene_id_known'], dropna=dropna).count().reset_index()
    temp2.rename({'gene_id_known': 'gene_id', 'Start': 'n_ss'}, axis=1, inplace=True)
    temp2.head()

    # number each gene as 1' 2' or 3'
    temp2 = temp2.sort_values(by='n_ss', ascending=False)
    temp2['rank'] = temp2.sort_values(by='n_ss', ascending=False).groupby(['Name']).cumcount()+1
    temp2.head()

    # number of genes w/ intersecting sss / ic
    temp3 = temp2[['Name', 'gene_id']].groupby('Name', dropna=dropna).nunique(dropna=dropna).reset_index()
    temp3.rename({'gene_id': 'n_genes'}, axis=1, inplace=True)
    temp3.head()

    temp5 = temp.merge(temp2, on='Name')
    temp5 = temp5.merge(temp3, on='Name')
    temp5['perc'] = (temp5['n_ss']/temp5['n_total_ss'])*100
    return f_ss_ic_df, temp5, fusion_df

def fix_talon_fusion_transcripts(talon_filt_ab,
                                 talon_gtf,
                                 ref_ics,
                                 ref_gtf,
                                 wc,
                                 ofile):
    """
    Fix fusion gene assignments based on the ss concordance between
    annotated (v29 or vM21) transcripts.
    * If a "fusion" transcript matches 2+ genes, keep it as fusion and
      label which genes it's a fusion of
    * If a "fusion" transcripts matches 1 gene, assign it to that
      gene instead of its novel gene
    * If a "fusion" transcript is entirely novel, remove it. These
      are 100% of the time monoexonic and this function will throw
      an error if that's not the case

    Additionally remove transcripts that are labelled as intergenic
    and are monoexonic.
    """
    if wc['species'] == 'human':
        ver = 'v40_cerberus'
    elif wc['species'] == 'mouse':
        ver = 'vM25_cerberus'
    temp, df, fusion_df = get_fusion_ss_table(talon_filt_ab, talon_gtf, ref_ics, 'v40_cerberus', include_novel=False)

    # 0-gene ss intersection transcripts
    novel_tids = list(set(fusion_df.Name.tolist())-set(temp.Name.tolist()))

    gtf_df = pr.read_gtf(talon_gtf, duplicate_attr=True, rename_attr=True)
    gtf_df = gtf_df.df
    gtf_df['gid_stable'] = cerberus.get_stable_gid(gtf_df, 'gene_id')

    ref_gtf_df = pr.read_gtf(ref_gtf, duplicate_attr=True, rename_attr=True)
    ref_gtf_df = ref_gtf_df.df
    ref_gtf_df['gid_stable'] = cerberus.get_stable_gid(ref_gtf_df, 'gene_id')

    # all ics w/ 2+ genes intersected -- fusion, keep their novelties
    # add an entry to reflect that they're fusions between the genes that they intersected with
    temp = df.loc[df.n_genes>=2]
    temp = temp[['Name', 'gene_id']].groupby('Name').agg(','.join).reset_index()
    temp.rename({'Name': 'transcript_id', 'gene_id': 'fusion_genes'}, axis=1, inplace=True)
    gtf_df = gtf_df.merge(temp, how='left', on='transcript_id')

    # all ics w/ just 1 gene intersected -- not fusion, update gene ids to intersected gene ids
    temp = df.loc[df.n_genes == 1]
    gene_cols = ['Source', 'gene_id', 'gene_name', 'gene_status', 'gene_type', 'talon_gene', 'havana_gene', 'level',
                 'antisense_gene', 'gene_antisense_to_IDs', 'intergenic_novel', 'fusion_novel']
    talon_gene_cols = ['gene_status', 'talon_gene', 'antisense_gene', 'gene_antisense_to_IDs', 'intergenic_novel', 'fusion_novel']
    for ind, entry in temp.iterrows():
        t = entry.Name
        inds = gtf_df.loc[gtf_df.transcript_id==t].index
        g = entry.gene_id
        dummy_gene_entry = gtf_df.loc[(gtf_df.gid_stable==g)&(gtf_df.Feature=='gene')]
        if len(dummy_gene_entry.index) == 0:
            # have to pull from reference gtf instead and add corresponding gene entry
            dummy_gene_entry = ref_gtf_df.loc[(ref_gtf_df.gid_stable==g)&(ref_gtf_df.Feature=='gene')]
            dummy_gene_entry[talon_gene_cols] = np.nan
            gtf_df = pd.concat([gtf_df, dummy_gene_entry], axis=0)
        assert len(dummy_gene_entry.index) == 1
        for c in gene_cols:
            gtf_df.loc[inds, c] = dummy_gene_entry[c].values[0]

    # all ics w/o any genes intersected -- change to intergenic genes? (need to double check)
    temp2 = fusion_df.loc[(fusion_df.Name.isin(novel_tids))&(fusion_df.Coordinates!='-')]
    assert len(temp2) == 0
    novel_gids = fusion_df.loc[(fusion_df.Name.isin(novel_tids)), 'gene_id']
    gtf_df = gtf_df.loc[~gtf_df.gid_stable.isin(novel_gids)]
    # jk delete these cause they're all genomic

    # get rid of monoexonic novel transcripts
    talon_df = pd.read_csv(talon_filt_ab, sep='\t')
    gene_novs = ['Antisense', 'Intergenic', 'Fusion']
    monoex_int_tids = talon_df.loc[(talon_df.n_exons==1)&(talon_df.gene_novelty.isin(gene_novs)), 'annot_transcript_id']
    gtf_df = gtf_df.loc[~gtf_df.transcript_id.isin(monoex_int_tids)]

    # drop stable gid, sort, update ends, and dump
    gtf_df.drop('gid_stable', axis=1, inplace=True)
    gtf_df = cerberus.sort_gtf(gtf_df)
    # mainly ripped out of cerberus code
    gtf_temp = gtf_df.copy(deep=True)
    for mode in ['tss', 'tes']:
        fwd, rev = cerberus.get_stranded_gtf_dfs(gtf_temp)
        df = pd.DataFrame()
        for strand, temp in zip(['+', '-'], [fwd, rev]):

            # fix gene boundaries
            temp = cerberus.update_gene_ends(temp, mode, strand)
            df = pd.concat([df, temp], ignore_index=True)

        gtf_temp = df.copy(deep=True)
    pr.PyRanges(gtf_temp).to_gtf(ofile)
