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

from .utils import *

def get_lr_bulk_sample_colors():
    c_dict, order = get_tissue_age_colors()

    # c2c12
    c_dict['c2c12_myoblast'] = '#ca79a7'
    c_dict['c2c12_myotube'] = '#009c73'
    order += ['c2c12_myoblast', 'c2c12_myotube']

    # forelimb
    c_dict['forelimb_e11'] = '#99ebec'
    c_dict['forelimb_e13'] = '#01ced0'
    order += ['forelimb_e11', 'forelimb_e13']


    # adrenal, hc, ctx
    # for t in ['adrenal', 'hippocampus', 'cortex']:
    #     c_dict[t] = get_tissue_colors()[0][t]
    # manually grabbed shades from shade picker to be one darker than
    # the last timecourse pt
    c_dict['adrenal_gland'] = '#8e361b'
    c_dict['hippocampus'] = '#9a3c4f'
    c_dict['cortex'] = '#634273'

    order += ['adrenal_gland', 'hippocampus', 'cortex']


    # f1219
    c_dict['f1219'] = '#4340bc'
    order += ['f1219']

    return c_dict, order

def get_tissue_age_colors():
    c_dict, order = get_tissue_colors()

    # get different color for each age / tissue
    new_c_dict = {}
    min_color = '#FFFFFF'
    ages = ['4d', '10d', '14d', '25d', '36d', '2mo', '18-20mo']
    order = []
    for tissue, max_color in c_dict.items():
        cmap = mpl.colors.LinearSegmentedColormap.from_list(tissue, [min_color, max_color], N=8)
        for i, age in enumerate(ages):
            key = '{}_{}'.format(tissue, age)
            new_c_dict[key] = mpl.colors.to_hex(cmap(i+1))
            order.append(key)

    return new_c_dict, order

def get_tissue_colors(cats=None, rgb=False):
    d = os.path.dirname(__file__)
    fname = f'{d}/../figures/ref/mouse/tissue_colors.tsv'
    df = pd.read_csv(fname, sep='\t')
    c_dict = {}
    for ind, entry in df.iterrows():
        c_dict[entry.tissue] = entry.color
    order = ['adrenal', 'hippocampus', 'cortex', 'gastroc', 'heart']

    if cats:
        pop_list = []
        for key in keys:
            if key not in cats:
                pop_list.append(key)
        for p in pop_list:
            del c_dict[p]
        order = [o for o in order if o in cats]

    if rgb:
        for key, item in c_dict.items():
            item = item[1:]
            r,g,b = tuple(int(item[i:i+2], 16) for i in (0, 2, 4))
            c_dict[key] = (r,g,b)

    return c_dict, order

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

def plot_human_sample_legend(swan_file, ofile):
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    c_dict, order = get_biosample_colors(species='human')

    df = swan.read(swan_file).adata.obs.copy(deep=True)
    df = df[['sample', 'sample_display', 'tissue_or_cell_line']]
    df = df.drop_duplicates()
    df = df.sort_values(by=['tissue_or_cell_line', 'sample_display'], ascending=True)
    order = df.sample_display.tolist()
    df['number'] = [i for i in range(len(df.index))]
    c_dict_2 = {}
    for key, item in c_dict.items():
        key2 = df.loc[df['sample'] == key, 'sample_display'].values[0]
        c_dict_2[key2] = item
    df['sample_display'] = df['sample_display'].astype('category')
    df['sample_display'].cat.categories = order
    import matplotlib.patches as patches
    samples = []
    for s in df.sample_display.cat.categories:
        c = c_dict_2[s]
        samples.append(patches.Rectangle((0,0),1,1,facecolor=c))
    ax = sns.scatterplot()
    plt.legend(samples, df.sample_display.cat.categories)

    ax.set_xticklabels('')
    ax.set_yticklabels('')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False)
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False)

    plt.savefig(ofile, dpi=700, bbox_inches='tight')

def plot_mouse_sample_legend(swan_file, ofile):
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    df = swan.read(swan_file).adata.obs.copy(deep=True)

    _, order = get_lr_bulk_sample_colors()
    c_dict, _ = get_biosample_colors('mouse')
    df['sample'] = df['sample'].astype('category')
    sample_order = df['sample'].cat.categories.tolist()
    sample_colors = [c_dict[s] for s in sample_order]
    order = [o for o in order if o in sample_order]

    df = df[['sample', 'sample_display']].drop_duplicates()
    c_dict_2 = {}
    for key, item in c_dict.items():
        try:
            key2 = df.loc[df['sample'] == key, 'sample_display'].values[0]
            c_dict_2[key2] = item
        except:
            pass
    # order = adata.obs['sample'].cat.categories
    display_order = [df.loc[df['sample']==s, 'sample_display'].values[0] for s in order]
    df['sample_display'] = df['sample_display'].astype('category')
    import pdb; pdb.set_trace()
    df['sample_display'].cat.categories = display_order
    import matplotlib.patches as patches
    samples = []
    for s in display_order:
        c = c_dict_2[s]
        samples.append(patches.Rectangle((0,0),1,1,facecolor=c))
    ax = sns.scatterplot()
    plt.legend(samples, df.sample_display.cat.categories)

    ax.set_xticklabels('')
    ax.set_yticklabels('')
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False)
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.savefig(ofile, dpi=700, bbox_inches='tight')

def human_undet_gene_go(ab, gene_subset, min_tpm, ofile):

    # get detected genes
    df = pd.read_csv(ab, sep='\t')
    df, gids = get_tpm_table(df,
                       how='gene',
                       min_tpm=min_tpm,
                       gene_subset=gene_subset)
    df['detected'] = True

    gene_df, b_counts, b_cat_counts = get_gtf_info(how='gene')
    gene_df.reset_index(inplace=True)
    gene_df['gid'] = cerberus.get_stable_gid(gene_df, 'gid')

    df = df.merge(gene_df, how='outer', left_index=True, right_on='gid')
    df.detected = df.detected.fillna(False)

    print(len(df.index))
    df = df.loc[df.biotype_category == 'protein_coding']
    print(len(df.index))

    dbs = ['GO_Biological_Process_2021',
           'GO_Cellular_Component_2021',
           'GO_Molecular_Function_2021',
           'KEGG_2021_Human',
           'KEGG_2019_Human']
    bm = gp.parser.Biomart()
    datasets = bm.get_datasets(mart='ENSEMBL_MART_ENSEMBL')
    datasets.loc[datasets.Description.str.contains('Human')]

    gids = df.loc[~df.detected, 'gid'].str.rsplit('.', n=1, expand=True)[0].to_frame()
    print(len(gids))
    gids = gids.squeeze().str.strip().tolist()
    gids = bm.query(dataset='hsapiens_gene_ensembl',
               attributes=['ensembl_gene_id', 'external_gene_name'],
               filters={'ensembl_gene_id': gids})
    gids = gids.loc[~gids.external_gene_name.isna()]
    gnames = gids.external_gene_name.squeeze().str.strip().tolist()
    go = gp.enrichr(gene_list=gnames,
                    gene_sets=dbs,
                    organism='Human',
                    description='undet_genes',
                    outdir='undet_genes_GO',
                    cutoff=0.5)

    def rm_go_number(df):
        df['term'] = df['Term'].str.split('\(GO', expand=True)[0]
        df['term'] = df.term.str.capitalize()
        return df

    # GO on undetected protein-coding genes
    df = pd.read_csv('undet_genes_GO/GO_Biological_Process_2021.Human.enrichr.reports.txt', sep='\t')
    n = 3
    df = df.head(n)
    df = rm_go_number(df)
    color = get_talon_nov_colors()[0]['Known']


    sns.set_context('paper', font_scale=2.2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    plt.figure(figsize=(3,4))
    ax = sns.barplot(data=df, x='Combined Score', y='term', color=color, saturation=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set(title='GO Molecular Component')

    xlabel = 'Enrichr Combined Score'
    ylabel = 'Undetected genes'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)
    plt.savefig(ofile, dpi=500, bbox_inches='tight')


def plot_perc_mane_det_by_len(ab,
                              filt_ab,
                              t_metadata,
                              min_gene_tpm,
                              min_tpm,
                              obs_col,
                              max_t_len,
                              figsize=(6,4),
                              fname='figures/mane_det_by_len.pdf'):
    # get detected genes
    g_df = pd.read_csv(ab, sep='\t')
    g_df = get_det_table(g_df,
                     how='gene',
                     min_tpm=min_gene_tpm,
                     gene_subset='polya',
                     groupby='library')
    gids = g_df.columns.tolist()

    # get all the mane transcripts
    # that are from genes we detect > 10 TPM
    t_df = pd.read_csv(t_metadata, sep='\t')
    t_df = t_df.loc[t_df.MANE_Select==True]
    t_df['gid_stable'] = cerberus.get_stable_gid(t_df, 'gid')
    t_df = t_df.loc[t_df.gid_stable.isin(gids)]

    # get all detected transcripts
    df = pd.read_csv(filt_ab, sep='\t')
    df = get_det_table(df,
                   how='iso',
                   min_tpm=min_tpm,
                   gene_subset='polya',
                   groupby=obs_col)
    tids = df.columns
    det_df = pd.DataFrame(data=tids, columns=['tid'])
    det_df['detected'] = True

    # merge detected transcripts in w/ mane transcripts
    t_df = t_df.merge(det_df, on='tid', how='left')
    t_df = t_df[['tid', 'gname', 't_len', 'detected']]
    t_df.detected.fillna(False, inplace=True)

    # only get shorter bois and bin by kb
    t_df = t_df.loc[t_df.t_len<=max_t_len]
    bins = [i for i in range(0, t_df.t_len.max()+1000, 1000)]
    t_df['len_bin'] = pd.cut(t_df.t_len, bins)

    # total transcripts / bin
    total_df = t_df[['tid', 'len_bin']].groupby('len_bin').count().reset_index()
    total_df.rename({'tid': 'n_total'}, axis=1, inplace=True)

    # total transcripts / bin by det status
    det_df = t_df[['tid', 'detected', 'len_bin']].groupby(['detected', 'len_bin']).count().reset_index()
    det_df = det_df.merge(total_df, on='len_bin', how='left')

    # calculate % detected and limit to only the detected transcripts
    det_df['perc'] = (det_df.tid/det_df.n_total)*100
    det_df.perc.fillna(0, inplace=True)
    det_df = det_df.loc[det_df.detected]

    # make nicer bin names
    hr_bin_dict = {}
    for b in det_df.len_bin.unique():
        bmin = str(int(int(str(b)[1:].split(',')[0])/1000))
        bmax = str(int(int(str(b).split(', ')[1][0:-1])/1000))
        hr_bin_dict[b] = f'{bmin}-{bmax}'
    det_df['hr_len_bin'] = det_df['len_bin'].map(hr_bin_dict)

    # make the beautiful plot
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    height = figsize[1]
    width = figsize[0]
    aspect = width/height

    color = get_talon_nov_colors()[0]['Known']

    ax = sns.catplot(det_df, x='hr_len_bin', y='perc', kind='bar',
                     color=color, saturation=1,
                     height=height, aspect=aspect)
    xlabel = 'Transcript length (kb)'
    ylabel = '% of detected GENCODE \n v40 MANE transcripts'
    ylim = (0,100)
    ax.set(xlabel=xlabel, ylabel=ylabel, ylim=ylim)
    ax.tick_params(axis="x", rotation=45)

    def add_perc_2(ax):
        ylim = ax.get_ylim()[1]
        n_cats = len(ax.patches)
        for p in ax.patches:
            perc = p.get_height()
            tot = det_df.loc[det_df.perc==perc, 'tid'].values[0]
            x = p.get_x() + p.get_width() / 2 - (0.02)*n_cats
            y = p.get_y() + p.get_height() + ylim*0.04
            ax.annotate(tot, (x, y), size = 18, rotation=90)

    a = ax.axes[0,0]
    add_perc_2(a)

    plt.savefig(fname, dpi=500, bbox_inches='tight')

    return det_df
