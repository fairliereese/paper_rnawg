import cerberus
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import upsetplot
from scipy import stats
from matplotlib.ticker import ScalarFormatter
import swan_vis as swan
import ternary
from sklearn import preprocessing
import pylab as pl
import matplotlib.ticker as tck
from collections import defaultdict
import plotly.graph_objects as go
import math
import plotly.io as pio
from pandarallel import pandarallel
import gseapy as gp
import os
import sys

p = os.getcwd()
sys.path.append(p)

from utils import *

def get_tissue_cell_line_colors(cats=None):
    tissue = '#e39f24'
    cell_line = '#7680e8'
    c_dict = {'cell_line': cell_line,
              'tissue': tissue}
    order = ['cell_line', 'tissue']
    c_dict, order = rm_color_cats(c_dict, order, cats)
    return c_dict, order

def get_ad_colors():
    c_dict = {'healthy': '#bb8f8f',
              'AD': '#b5bd61',
              np.nan: '#000000'}
    order = ('healthy', 'AD', np.nan)
    return c_dict, order

def get_not_det_color():
    return '#E5ECF6'

def get_shade_colors(color, order):
    c_dict = {}
    min_color = '#FFFFFF'
    cmap = mpl.colors.LinearSegmentedColormap.from_list('temp', [color, min_color], N=len(order)+1)
    for i, cat in enumerate(order):
        c_dict[cat] = mpl.colors.to_hex(cmap(i))

    return c_dict, order

def get_biosample_colors(species='human'):
    """
    Get colors for each biosample
    """
    d = os.path.dirname(__file__)
    fname = f'{d}/data/{species}/lr/lr_{species}_library_data_summary.tsv'
    df = pd.read_csv(fname, sep='\t')
    df = df[['sample', 'sample_color_hex_code']].drop_duplicates()

    c_dict = {}
    for ind, entry in df.iterrows():
        c_dict[entry['sample']] = entry.sample_color_hex_code
    order = df['sample'].tolist()

    return c_dict, order

def rm_color_cats(c_dict, order, cats):
    if cats:
        keys = c_dict.keys()
        pop_list = []
        for key in keys:
            if key not in cats:
                pop_list.append(key)
        for p in pop_list:
            del c_dict[p]
        order = [o for o in order if o in cats]
    return c_dict, order

def get_ic_nov_colors(cats=None):
    c_dict = {'Known': '#009E73',
              'ISM': '#0072B2',
              'NIC': '#D55E00',
              'NNC': '#E69F00',
              'Antisense': '#000000',
              'Intergenic': '#CC79A7',
              'Unspliced': '#F0E442'}
    order = ['Known', 'ISM', 'NIC', 'NNC', 'Antisense', 'Intergenic', 'Unspliced']

    c_dict, order = rm_color_cats(c_dict, order, cats)
    return c_dict, order

def plot_gene_det_by_biotype_tpm(df,
                                 how,
                                 ver,
                                 opref='figures/',
                                 **kwargs):
    if how == 'sr':
        df, inds = get_tpm_table(df,
                    how='sr',
                    gene_subset='polya',
                    min_tpm=0,
                    **kwargs)
    else:
        df, inds = get_tpm_table(df,
                    how='gene',
                    gene_subset='polya',
                    min_tpm=0,
                    **kwargs)

    gene_df, b_counts, b_cat_counts = get_gtf_info(how='gene', ver=ver)

    polya_biotypes = ['protein_coding', 'pseudogene', 'lncRNA']
    polya_genes = gene_df.loc[gene_df.biotype_category.isin(polya_biotypes), 'gid'].tolist()
    n_polya = len(polya_genes)
    n_det_polya = len(df.index)

    print('Detected {} / {} ({:.3}%) annotated polyA genes'.format(n_det_polya, n_polya, (n_det_polya/n_polya)*100))

    tpm_df = df.copy(deep=True)
    tpm_dfs = []
    tpm_dfs.append(tpm_df)
    tpm_dfs.append(tpm_df.loc[(tpm_df >= 1).any(axis=1)])
    tpm_dfs.append(tpm_df.loc[(tpm_df >= 100).any(axis=1)])

    det_df = pd.DataFrame()
    for df, tpm in zip(tpm_dfs, [0,1,100]):
        gene_df, _, _ = get_gtf_info(how='gene', subset='polya', ver=ver, add_stable_gid=True)
        gene_df = gene_df[['gid_stable', 'gname', 'biotype_category']]

        df.reset_index(inplace=True)
        df.rename({'index': 'gid_stable'}, axis=1, inplace=True)
        df = df.merge(gene_df, how='left', on='gid_stable')

        df = df[['gid_stable', 'biotype_category']].groupby('biotype_category').count()
        df.rename({'gid_stable':'obs_counts'}, axis=1, inplace=True)

        gene_df = gene_df[['gid_stable', 'biotype_category']].groupby('biotype_category').count()
        gene_df.rename({'gid_stable':'annot_counts'}, axis=1, inplace=True)
        df = df.merge(gene_df, how='left', left_index=True, right_index=True)

        df['perc'] = (df.obs_counts/df.annot_counts)*100
        df = df.sort_values(by='perc', ascending=False)
        df['tpm_thresh'] = tpm
        det_df = pd.concat([det_df, df])

    det_df = det_df.reset_index()
    det_df.head()

    sns.set_context('paper', font_scale=2)
    plt.figure(figsize=(6,6))
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    ic_colors, order = get_ic_nov_colors()
    gray = get_not_det_color()
    c = ic_colors['Known']
    cats = [100,1,0]
    c_dict, order = get_shade_colors(c, cats)
    order.reverse()
    biotypes = ['protein_coding', 'lncRNA', 'pseudogene']
    b_dict = {'protein_coding': 'Protein\ncoding',
              'lncRNA': 'lncRNA',
              'pseudogene': 'Pseudogene'}

    # https://matplotlib.org/2.0.2/examples/api/barchart_demo.html
    def add_n(rects, label):
        ax = plt.gca()
        for rect in rects:
            # height = rect.get_height()
            x = rect.get_y()+rect.get_height()/2.5
            y = rect.get_width()*1.1
            ax.text(y,x,
                    '{:,}'.format(label),
                    ha='center', va='bottom', size=16)

    def add_n_2(rects, label):
        ax = plt.gca()
        for rect in rects:
            # height = rect.get_height()
            x = rect.get_y()+rect.get_height()*1.2
            y = rect.get_width()/2
            ax.text(y,x,
                    '{:,}'.format(label),
                    ha='center', va='bottom', size=16)

    for b in biotypes:
        x = b_dict[b]
        y = 0
        rects = plt.barh(x, [100], color=gray, height=0.5)
        # add total number of genes
        n = det_df.loc[(det_df.biotype_category == b)&(det_df.tpm_thresh==1), 'annot_counts'].tolist()[0]
        add_n(rects, n)

        for c in order:
            curr_y = det_df.loc[(det_df.biotype_category == b)&(det_df.tpm_thresh==c), 'perc'].tolist()[0]
            rects = plt.barh(x, [curr_y], color=c_dict[c], height=0.5)
            if c == 1:
                n = det_df.loc[(det_df.biotype_category == b)&(det_df.tpm_thresh==1), 'obs_counts'].tolist()[0]
                add_n_2(rects, n)
                print(b)
                print(curr_y)
                print(det_df.loc[(det_df.biotype_category == b)&(det_df.tpm_thresh==c), ['obs_counts', 'annot_counts']])
                print()
            y = y+curr_y



    leg_labels = ['Not Detected', 'Detected', 'Detected >= 1 TPM', 'Detected >= 100 TPM']
    plt.legend(leg_labels, bbox_to_anchor=(.6, 1.05))
    ax = plt.gca()
    leg = ax.get_legend()

    # plt.yticks(rotation=90)

    # plt.ylabel('Biotype')
    plt.xlabel('% of GENCODE v40 genes')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    fname = f'{opref}gene_det_by_biotype.png'
    plt.savefig(fname, dpi=500, bbox_inches='tight')
    fname = f'{opref}gene_det_by_biotype.pdf'
    plt.savefig(fname, dpi=500, bbox_inches='tight')

def plot_species_sector_gene_counts(m_counts, h_counts):
    temp = pd.DataFrame()
    for source in ['GENCODE', 'obs']:
        for species, counts in zip(['mouse', 'human'],[m_counts, h_counts]):
            df = assign_gisx_sector(counts)
            df = df.loc[df.source == source]
            df = df[['gid', 'source', 'sector']].groupby(['source', 'sector']).count().reset_index()
            df.rename({'gid': 'n_genes'}, axis=1, inplace=True)
            df['total_genes'] = df.n_genes.sum()
            df['species'] = species
            temp = pd.concat([temp, df])
    temp['perc'] = (temp.n_genes/temp.total_genes)*100

    y = '% annotated / observed genes'
    temp.rename({'perc': y}, axis=1, inplace=True)
    c_dict, order = get_sector_colors(['tss', 'splicing', 'tes'])
    temp = temp.loc[temp.sector != 'simple']

    # plot both together
    sns.set_context('paper', font_scale=1.8)
    ax = sns.catplot(data=temp, x='source',
                y=y, col='species',
                hue='sector', kind='bar',
                hue_order=order,
                palette=c_dict, saturation=1)

    def add_perc_2(ax):
        ylim = ax.get_ylim()[1]
        n_cats = len(ax.patches)
        for p in ax.patches:
            percentage = '{:.1f}%'.format(p.get_height())
    #         x = p.get_x() + p.get_width() / 2 - 0.45
            x = p.get_x() + p.get_width() / 2 - (0.015)*n_cats
            y = p.get_y() + p.get_height() + ylim*0.00625
            ax.annotate(percentage, (x, y), size = 12)

    a = ax.axes[0,0]
    add_perc_2(a)
    a = ax.axes[0,1]
    add_perc_2(a)

    return temp

def plot_biosamp_det(df,
                     opref='figures/',
                     figsize=(6,6),
                     **kwargs):

    """
    Plot a hist of the number of tissues, cell lines, or datasets each gene or
    isoform occurs in.

    Parameters:
        df (pandas DataFrame): TALON abundance
        opref (str): Output prefix to save figure

    Returns:
        df (pandas DataFrame): DataFrame detailing how many samples
            each gene or isoform was seen in
    """
    if 'ic_nov' in kwargs:
        nov = kwargs['ic_nov'][0]
    if 'how' in kwargs:
        how = kwargs['how']
    if 'groupby' in kwargs:
        groupby = kwargs['groupby']
    if 'sample' in kwargs:
        sample = kwargs['sample']
    if 'nov' in kwargs:
        nov = kwargs['nov'][0]
    else:
        nov = 'Known'

    df = get_det_table(df, **kwargs)

    # finally, calculate the number of biosamples / libraries these
    # genes or transcripts are expressed >= min TPM
    df = df.transpose()
    df['n_samples'] = df.astype(int).sum(axis=1)

    # and make a beautiful plot
    plt.figure(figsize=figsize)
    height = figsize[1]
    width = figsize[0]
    aspect = width/height
    # print(height)
    # print(aspect)
    # https://stackoverflow.com/questions/65415646/how-to-change-the-figure-size-of-a-displot
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    c_dict, order = get_ic_nov_colors()
    if nov:
        color = c_dict[nov]
    else:
        color = c_dict['Known']

    # color by triplet feature type
    if how in ['tss', 'tes', 'ic']:
        c_dict, order = get_sector_colors()
        if how == 'ic':
            color = c_dict['splicing']
        else:
            color = c_dict[how]

    ax = sns.displot(data=df, x='n_samples', kind='hist',
                 color=color, binwidth=1, linewidth=0,
                 alpha=1,height=height, aspect=aspect)

    # titles
    if how == 'gene' or how == 'sr':
        ylabel = '# known genes'
    elif how == 'iso':
        if nov == 'Known':
            nov = nov.lower()
        ylabel = '# transcripts w/ {} ICs'.format(nov)
    elif how in ['tss', 'ic', 'tes']:
        ylabel = '# {}s'.format(how.upper())

    if groupby == 'sample':
        xlabel = '# samples'
    elif groupby == 'biosample':
        xlabel = '# ENCODE biosamples'
    elif groupby == 'tissue':
        xlabel = '# tissues'
    elif groupby == 'cell_line':
        xlabel = '# celltypes'
    elif groupby == 'library':
        if sample == 'cell_line':
            xlabel = '# cell line libraries'
        elif sample == 'tissue':
            xlabel = '# tissue libraries'
        else:
            xlabel = '# libraries'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    sample = groupby

    if groupby == 'sample':
        fname = '{}{}_{}_{}_detection.png'.format(opref, sample, nov, how)
    else:
        fname = '{}{}_{}_{}_library_detection.png'.format(opref, sample, nov, how)

    plt.savefig(fname, dpi=500, bbox_inches='tight')

    if groupby == sample:
        fname = '{}{}_{}_{}_detection.pdf'.format(opref, sample, nov, how)
    else:
        fname = '{}{}_{}_{}_library_detection.pdf'.format(opref, sample, nov, how)

    plt.savefig(fname, dpi=500, bbox_inches='tight')

    return df
