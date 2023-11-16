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

# from .utils import *
from utils import *

def get_sector_colors(cats=None):
    tss = '#56B4E9'
    tes = '#E69F00'
    splicing = '#CC79A7'
    simple = '#000000'
    c_dict = {'tss': tss,
              'splicing': splicing,
              'tes': tes,
              'simple': simple,
              'mixed': '#b7b7b7'}
    order = ['tss', 'splicing', 'tes', 'mixed', 'simple']

    c_dict, order = rm_color_cats(c_dict, order, cats)
    return c_dict, order

def get_talon_nov_colors(cats=None):
    c_dict = {'Known': '#009E73',
              'ISM': '#0072B2',
              'ISM_rescue': '#0072B2',
              'NIC': '#D55E00',
              'NNC': '#E69F00',
              'Antisense': '#000000',
              'Intergenic': '#CC79A7',
              'Genomic': '#F0E442'}
    order = ['Known', 'ISM', 'ISM_rescue', 'NIC', 'NNC', 'Antisense', 'Intergenic', 'Genomic']

    c_dict, order = rm_color_cats(c_dict, order, cats)
    return c_dict, order

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


    gids = df.loc[~df.detected, 'gid'].str.rsplit('.', n=1, expand=True)[0].to_frame()
    df = pd.DataFrame()
    df['gid'] = gids
    gene_df, _, _ = get_gtf_info(how='gene', ver='v40_cerberus', add_stable_gid=True)
    df = df.merge(gene_df[['gid_stable', 'gname']], how='left', left_on='gid', right_on='gid_stable')
    gnames = df.gname.tolist()

    go = gp.enrichr(gene_list=gnames,
                    gene_sets=dbs,
                    organism='Human',
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
    # ax.set(title='GO Molecular Component')

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
                     groupby='library',
                     species='human')
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

def plot_mouse_cell_line_tissue_read_len_v_ref(read_annot,
                                         meta_file,
                                         gene_subset,
                                         ref_fname,
                                         xlim,
                                         ofile):

    df = pd.read_csv(read_annot, usecols=[1,8], sep='\t')
    meta = pd.read_csv(meta_file, sep='\t')
    cell_lines = meta.loc[meta.tissue_or_cell_line == 'cell_line', 'dataset'].tolist()
    tissues = meta.loc[meta.tissue_or_cell_line == 'tissue', 'dataset'].tolist()
    df['source'] = False
    df.loc[df.dataset.isin(cell_lines), 'source'] = 'Reads from cell lines'
    df.loc[df.dataset.isin(tissues), 'source'] = 'Reads from tissues'

    t_df, _, _ = get_gtf_info('iso',
                              fname=ref_fname,
                              subset=gene_subset)

    df = df[['read_length', 'source']]
    df.rename({'read_length': 'length'}, axis=1, inplace=True)
    # df['source'] = 'Reads'
    t_df = t_df[['t_len']]
    t_df.rename({'t_len':'length'}, axis=1, inplace=True)
    t_df['source'] = 'GENCODE vM25 transcripts'
    df = pd.concat([df, t_df])

    sns.set_context('paper', font_scale=2)
    plt.figure(figsize=(3,3))
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    temp_c_dict, order = get_ic_nov_colors()
    c_dict = dict()
    c_dict['GENCODE vM25 transcripts'] = temp_c_dict['Known']
    temp_c_dict, order = get_tissue_cell_line_colors()
    c_dict['Reads from tissues'] = temp_c_dict['tissue']
    c_dict['Reads from cell lines'] = temp_c_dict['cell_line']
    order = ['Reads from cell lines',
             'Reads from tissues',
             'GENCODE vM25 transcripts']

    ax = sns.displot(data=df, x='length', kind='kde',
                         linewidth=3, common_norm=False, hue='source',
                         palette=c_dict, hue_order=order)
    xlabel = 'Read / transcript length'
    ylabel = 'Density'

    if xlim:
        _ = ax.set(xlabel=xlabel, ylabel=ylabel, xlim=(0,xlim))
    else:
        _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    plt.savefig(ofile, dpi=500, bbox_inches='tight')

def plot_cell_line_tissue_read_len_v_ref(df,
                                         lib_meta,
                                         gene_subset,
                                         ref_fname,
                                         xlim,
                                         ofile):

    # cell_lines = get_sample_datasets('human', 'cell_line')
    # tissues = get_sample_datasets('human', 'tissue')
    df['source'] = False
    # df.loc[df.dataset.isin(cell_lines), 'source'] = 'Reads from cell lines'
    # df.loc[df.dataset.isin(tissues), 'source'] = 'Reads from tissues'
    meta = pd.read_csv(lib_meta, sep='\t')
    meta = meta[['dataset', 'tissue_or_cell_line']]
    df = df.merge(meta, how='left', on='dataset')
    df.loc[df.tissue_or_cell_line=='tissue', 'source'] = 'Reads from tissues'
    df.loc[df.tissue_or_cell_line=='cell_line', 'source'] = 'Reads from cell lines'

    t_df, _, _ = get_gtf_info('iso',
                          fname=ref_fname,
                          subset=gene_subset)

    df = df[['read_length', 'source']]
    df.rename({'read_length': 'length'}, axis=1, inplace=True)
    t_df = t_df[['t_len']]
    t_df.rename({'t_len':'length'}, axis=1, inplace=True)
    t_df['source'] = 'GENCODE v40 transcripts'
    df = pd.concat([df, t_df])

    sns.set_context('paper', font_scale=2)
    plt.figure(figsize=(1.5, 2))
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    temp_c_dict, order = get_ic_nov_colors()
    c_dict = dict()
    c_dict['GENCODE v40 transcripts'] = temp_c_dict['Known']
    temp_c_dict, order = get_tissue_cell_line_colors()
    c_dict['Reads from tissues'] = '#e39f24'
    c_dict['Reads from cell lines'] = '#7680e8'
    order = ['Reads from cell lines',
             'Reads from tissues',
             'GENCODE v40 transcripts']

    ax = sns.displot(data=df, x='length', kind='kde',
                         linewidth=3, common_norm=False, hue='source',
                         palette=c_dict, hue_order=order)
    xlabel = 'Read / transcript length'
    ylabel = 'Density'

    if xlim:
        _ = ax.set(xlabel=xlabel, ylabel=ylabel, xlim=(0,xlim))
    else:
        _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    plt.savefig(ofile, dpi=500, bbox_inches='tight')

def plot_supported_feats_2(filt_ab,
                         h5,
                         feat,
                         ref_sources,
                         support_sources,
                         opref,
                         ax,
                         **kwargs):
    """
    Plot a bar plot showing which observed features are supported
    from external data sources in the cerberus annotation

    Parameters:
        filt_ab (str): Path fo filtered abundance file
        h5 (str): Path to cerberus annotation h5 object
        feat (str): {'tss', 'tes', 'ic'}
    """

    # get detected features
    df = pd.read_csv(filt_ab, sep='\t')
    df, ids = get_tpm_table(df, **kwargs)

    # get these features from cerberus
    ca = cerberus.read(h5)
    if feat == 'tss':
        ca_df = ca.tss
    elif feat == 'tes':
        ca_df = ca.tes
    elif feat == 'ic':
        ca_df = ca.ic
    print(len(ca_df.index))
    df = ca_df.loc[ca_df.Name.isin(ids)]
    print(len(df.index))


    # get T/F detection of each feat by each source
    df = upsetplot.from_memberships(df.source.str.split(','), data=df)
    df.reset_index(inplace=True)

    # which sources are observed, which are supported, and which are known
    sources = ca.get_sources(df)

    df['support'] = 'Novel'
    if support_sources:
        df.loc[df[support_sources].any(axis=1), 'support'] = 'Supported'
    df.loc[df[ref_sources].any(axis=1), 'support'] = 'Known'
    df = df[['Name', 'support']]
    df = df.groupby('support').count().reset_index()
    df.rename({'Name': 'counts'}, axis=1, inplace=True)
    print(df)

    # colors
    # pdb.set_trace()
    if feat == 'ic':
        c_key = 'splicing'
    else:
        c_key = feat
    temp_c_dict, order = get_sector_colors()
    c = temp_c_dict[c_key]
    order = ['Known', 'Supported', 'Novel']
    c_dict, order = get_shade_colors(c, order)


    # plotting
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    # plt.figure(figsize=(3,4))


    # ax = plt.gca()

    # print out some shtuff
    n = df.counts.sum()
    n_num = df.loc[df.support!='Known'].counts.sum()
    print('{:.2f}% ({}/{}) of {}s are novel'.format((n_num/n)*100, n_num, n, feat))

    n = df.loc[df.support!='Known'].counts.sum()
    n_num = df.loc[df.support == 'Supported'].counts.sum()
    print('{:.2f}% ({}/{}) of novel {}s are supported'.format((n_num/n)*100, n_num, n, feat))


    # Known
    c = 'Known'
    y = df.loc[df.support==c, 'counts'].values[0]
    ax.bar(c, y, color=c_dict[c])

    # Novel
    c = 'Novel'
    nov_list = ['Novel', 'Supported']
    y = df.loc[df.support.isin(nov_list), 'counts'].sum(axis=0)
    ax.bar(c, y, color=c_dict[c])

    # Supported novel
    c = 'Supported'
    y = df.loc[df.support==c, 'counts'].values[0]
    ax.bar('Novel', y, color=c_dict[c])


    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = ''
    ylabel = '# observed {}s'.format(feat.upper())

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.tick_params(axis="x", rotation=45)

    # plt.legend()

    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[0] = 'Known'
    labels[1] = 'Novel'
    ax.set_xticklabels(labels)

    # leg_list = ['Novel', 'Supported', 'Known']
    # # if feat == 'tss':
    # #     support_label = "5' end support"
    # # elif feat == 'tes':
    # #     support_label = "3' end support"
    # # elif feat == 'ic':
    # #     support_label = "IC support"
    # support_label = 'Supported'
    # leg_list.append(support_label)
    # leg_list.reverse()

    ax.legend(order)
    leg = ax.get_legend()
    c_dict, order = get_shade_colors('#000000', order)
    # pdb.set_trace()
    # # print(c_dict)
    # c_dict, order = get_shade_colors(c_dict['Novel'], leg_list)
    for i, nov in enumerate(order):
        leg.legendHandles[i].set_color(c_dict[nov])

    fname = '{}_{}_support.png'.format(opref, feat)
    plt.savefig(fname, dpi=500, bbox_inches='tight', layout='tight')

    fname = '{}_{}_support.pdf'.format(opref, feat)
    plt.savefig(fname, dpi=500, bbox_inches='tight', layout='tight')
    print()

def plot_n_feat_per_gene(h5,
                         feat,
                         max_ends=10,
                         show_pc=False,
                         subset=None):
    """
    Plot number of features per gene in a given source,
    in a given subset
    """

    feat_col = 'n_{}'.format(feat)

    # get these features from cerberus
    df = get_ca_table(h5, feat)

    # get subset features
    if subset:
        df = df.loc[df.Name.isin(subset)]

    #
    if show_pc:
        gene_df, _, _ = get_gtf_info(how='gene', ver='v40_cerberus')
        gene_df['gid'] = cerberus.get_stable_gid(gene_df, col='gid')
        gene_df = gene_df[['gid', 'biotype_category']]
        df = df.merge(gene_df, how='left', left_on='gene_id', right_on='gid')

    # count # feats / gene
    gb_cols = ['gene_id']
    if show_pc:
        gb_cols += ['biotype_category']
    keep_cols = gb_cols + ['Name']
    df = df[keep_cols]
    df = df.groupby(gb_cols).count().reset_index()
    df.rename({'Name': feat_col}, axis=1, inplace=True)

    # create counts df
    gb_cols = [feat_col]
    if show_pc:
        gb_cols += ['biotype_category']
    df = df.groupby(gb_cols).count().reset_index()

    # group all the entries over the max number
    df.rename({'gene_id': 'n_genes'}, axis=1, inplace=True)
    df.loc[df[feat_col] >= max_ends, feat_col] = '{}+'.format(max_ends)
    df = df.groupby(gb_cols).sum().reset_index()

    # pdb.set_trace()

    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    plt.figure(figsize=(3,4))

    c_dict, order = get_sector_colors()
    if feat == 'ic':
        c = c_dict['splicing']
    else:
        c = c_dict[feat]
    if show_pc:
        biotypes = ['protein_coding', 'lncRNA', 'pseudogene']
        b_dict = {'protein_coding': 'Protein coding',
                  'lncRNA': 'lncRNA',
                  'pseudogene': 'Pseudogene'}
        # biotypes.reverse()
        c_dict, order = get_shade_colors(c, biotypes)
        # order.reverse()

    df = df.pivot(index=feat_col, columns=['biotype_category'])
    df.columns = df.columns.droplevel(0)
    df.reset_index(inplace=True)

    df[feat_col] = df[feat_col].astype(str)
    x = df[feat_col].unique().tolist()

    # loop through biotypes
    bottom = [0 for i in range(len(x))]
    for b in order:
        y = df[b].tolist()
        plt.bar(x, y, color=c_dict[b], bottom=bottom)
        bottom = [b_coord+y_coord for b_coord, y_coord in zip(bottom, y)]

    plt.xlabel('# {}s / gene'.format(feat.upper()))
    plt.ylabel('# genes')
    sns.despine()

    leg_labels = [b_dict[o] for o in order]
    plt.legend(leg_labels, bbox_to_anchor=(.6, 1.05))
    ax = plt.gca()
    leg = ax.get_legend()
    shade_dict, _ = get_shade_colors('#000000', order)
    for i, o in enumerate(order):
        leg.legendHandles[i].set_color(shade_dict[o])

    fname = 'figures/{}_per_gene_support.png'.format(feat)
    plt.savefig(fname, dpi=500, bbox_inches='tight')

    fname = 'figures/{}_per_gene_support.pdf'.format(feat)
    plt.savefig(fname, dpi=500, bbox_inches='tight')

    return df



def plot_novel_supported_triplet_feats(filt_ab,
                                       h5,
                                       gene_subset,
                                       min_tpm,
                                       ofile):

    feats = ['tss', 'ic', 'tes']
    ref_sources = [['v29', 'v40'],
                   ['v29', 'v40'],
                   ['v29', 'v40']]
    support_sources = [['encode_cage', 'fantom_cage', 'encode_rampage', 'gtex'],
                       ['gtex'],
                       ['pas', 'polya_atlas', 'gtex']]

    plt.figure(figsize=(3,20))
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    fig, axes = plt.subplots(1, len(feats), figsize=(10,4))
    i = 0

    for feat, ref_source, support_source in zip(feats, ref_sources, support_sources):
        print(feat)
        ax = axes[i]
        print(f'Support sources:{support_source}')
        plot_supported_feats_2(filt_ab=filt_ab,
                             h5=h5,
                             feat=feat,
                             ref_sources=ref_source,
                             support_sources=support_source,
                             how=feat,
                             opref='figures/human',
                             gene_subset=gene_subset,
                             min_tpm=1,
                             ax=ax)
        i += 1

    plt.subplots_adjust(wspace=0.35)

    plt.savefig(ofile, dpi=500, layout='tight', bbox_inches='tight')

def plot_triplet_feats_per_gene(h5,
                                filt_ab,
                                gene_subset,
                                min_tpm,
                                opref='figures/'):
    dfs = dict()
    for feat in ['tss', 'ic', 'tes']:
        ids = get_det_feats(h5,
                            filt_ab,
                            feat,
                            how=feat,
                            gene_subset=gene_subset,
                            min_tpm=min_tpm)
        df = plot_n_feat_per_gene(h5,
                             feat=feat,
                             show_pc=True,
                             subset=ids)
        dfs[feat] = df
        fname = opref+feat+'.pdf'
        plt.savefig(fname, dpi=500, layout='tight', bbox_inches='tight')
    return dfs

def plot_transcript_det_by_biotype(filt_ab,
                                   h5,
                                   min_tpm,
                                   gene_subset,
                                   ic_nov,
                                   species,
                                   ofile):

    df = pd.read_csv(filt_ab, sep='\t')
    df, tids = get_tpm_table(df,
                   how='iso',
                   min_tpm=min_tpm,
                   gene_subset=gene_subset,
                   ic_nov=ic_nov,
                   species=species,
                   h5=h5)
    df['temp_gid'] = df.index.to_series().str.split('[',expand=True)[0]

    n = len(df.index)
    print('Detected {} transcripts >= 1 TPM w/ a known ic'.format(n))

    gene_df, _, _ = get_gtf_info(how='gene', ver='v40_cerberus', subset='polya', add_stable_gid=True)
    df = df.merge(gene_df, how='left',
                  left_on='temp_gid',
                  right_on='gid_stable')

    df = df[['biotype_category']]
    df.reset_index(inplace=True)
    df = df.groupby('biotype_category').count().reset_index()

    c = get_ic_nov_colors()[0]['Known']
    c_dict, order = get_shade_colors(c, ['Protein\ncoding', 'lncRNA', 'Pseudogene'])
    order.reverse()

    label_dict = {'protein_coding': 'Protein\ncoding',
                  'lncRNA': 'lncRNA',
                  'pseudogene': 'Pseudogene'}

    df['label_biotype_cat'] = df.biotype_category.map(label_dict)

    # plotting
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    # plt.figure(figsize=(4,6))
    plt.figure(figsize=(3,4))

    ax = sns.barplot(data=df, y='label_biotype_cat', x='index',
                     palette=c_dict, order=order,
                     saturation=1)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # ylabel = 'RNA biotype'
    ylabel = ''
    xlabel = '# transcripts w/ known ICs'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.tick_params(axis="y")


    # fname = f'{fig_dir}/iso_det_by_biotype.pdf'
    # plt.savefig(fname, dpi=500, bbox_inches='tight')

    plt.savefig(ofile, dpi=500, bbox_inches='tight')


def plot_transcripts_by_triplet_feat_novelty(filt_ab,
                                             h5,
                                             min_tpm,
                                             gene_subset,
                                             species,
                                             ofile):

    df = pd.read_csv(filt_ab, sep='\t')
    df, tids = get_tpm_table(df,
                   how='iso',
                   min_tpm=min_tpm,
                   gene_subset=gene_subset,
                   species=species,
                   h5=h5)

    df = df.reset_index()
    for feat in ['tss', 'ic', 'tes']:
        df = add_feat(df, col='annot_transcript_id', kind=feat)
        feat_df = get_ca_table(h5, feat)
        feat_df = feat_df[['Name', 'novelty']]
        feat_df.rename({'novelty': '{}_novelty'.format(feat),
                        'Name': feat}, axis=1, inplace=True)
        df = df.merge(feat_df, how='left', on=feat)

    df = df[['annot_transcript_id', 'tss_novelty', \
             'ic_novelty', 'tes_novelty']].groupby(['tss_novelty', \
             'ic_novelty', \
             'tes_novelty']).count()
    df = df.reset_index()

    # determine novelty categories
    df['ends_novelty'] = np.nan
    inds = df.loc[(df.tss_novelty == 'Known')&(df.tes_novelty == 'Known')].index.tolist()
    df.loc[inds, 'ends_novelty'] = 'Known'
    inds = df.loc[(df.tss_novelty == 'Known')&(df.tes_novelty == 'Novel')].index.tolist()
    df.loc[inds, 'ends_novelty'] = 'Novel_tes'
    inds = df.loc[(df.tss_novelty == 'Novel')&(df.tes_novelty == 'Known')].index.tolist()
    df.loc[inds, 'ends_novelty'] = 'Novel_tss'
    inds = df.loc[(df.tss_novelty == 'Novel')&(df.tes_novelty == 'Novel')].index.tolist()
    df.loc[inds, 'ends_novelty'] = 'Novel'

    df.rename({'annot_transcript_id': 'counts'}, axis=1, inplace=True)
    n = df.counts.sum()
    n_num = df.loc[(df.tss_novelty=='Known')&(df.ic_novelty=='Known')&(df.tes_novelty=='Known')].counts.sum()
    print('{:.2f}% ({}/{}) transcripts have all known triplet features'.format((n_num/n)*100, n_num, n))

    ic_colors, order = get_ic_nov_colors(df.ic_novelty.unique().tolist())

    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    for nov in order:
        c = ic_colors[nov]
        c_dict, order = get_shade_colors(c, ['Known','Novel_tes','Novel_tss','Novel'])
        x = [nov]
        y = 0
        for end_novelty in order:
            # print(end_novelty)
            curr_y = df.loc[(df.ic_novelty == nov)&(df.ends_novelty==end_novelty), 'counts'].tolist()[0]
            # print('curr_y: ', curr_y)
            # print('y: ', y)
            # print()
            plt.bar(x, [curr_y], bottom=y, color=c_dict[end_novelty])
            y = y+curr_y
    plt.legend(['Known TSS+TES', 'Novel TSS', 'Novel TES', 'Novel TSS+TES'])
    ax = plt.gca()
    leg = ax.get_legend()
    c_dict, order = get_shade_colors('#000000', ['Known','Novel_tes','Novel_tss','Novel'])
    for i, nov in enumerate(order):
        leg.legendHandles[i].set_color(c_dict[nov])
    # plt.xlabel('IC novelty')
    plt.ylabel('# transcripts')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis="x", rotation=45)

    plt.savefig(ofile, dpi=500, bbox_inches='tight')

    return df

def plot_n_isos_per_gene(df,
                       max_isos=10,
                       show_pc=False,
                       subset=None,
                    opref='figures/'):

    """
    Plot number of isoforms per gene in a
    given subset of isoforms
    """

    df['gid'] = cerberus.get_stable_gid(df, 'annot_gene_id')

    feat_col = 'n_iso'

    # get subset features
    if subset:
        df = df.loc[df.annot_transcript_id.isin(subset)]

    #
    if show_pc:
        gene_df, _, _ = get_gtf_info(how='gene', ver='v40_cerberus')
        gene_df['gid'] = cerberus.get_stable_gid(gene_df, col='gid')
        gene_df = gene_df[['gid', 'biotype_category']]
        df = df.merge(gene_df, how='left', on='gid')

    # count # feats / gene
    gb_cols = ['gid']
    if show_pc:
        gb_cols += ['biotype_category']
    keep_cols = gb_cols + ['annot_transcript_id']
    df = df[keep_cols]
    df = df.groupby(gb_cols).count().reset_index()
    df.rename({'annot_transcript_id': feat_col}, axis=1, inplace=True)

    # create counts df
    gb_cols = [feat_col]
    if show_pc:
        gb_cols += ['biotype_category']
    df = df.groupby(gb_cols).count().reset_index()

    # group all the entries over the max number
    df.rename({'gene_id': 'n_genes'}, axis=1, inplace=True)
    df.loc[df[feat_col] >= max_isos, feat_col] = '{}+'.format(max_isos)
    df = df.groupby(gb_cols).sum().reset_index()

    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    plt.figure(figsize=(6,4))

    c_dict, order = get_talon_nov_colors(cats=['Known'])
    c = c_dict['Known']
    if show_pc:
        biotypes = ['protein_coding', 'lncRNA', 'pseudogene']
        b_dict = {'protein_coding': 'Protein coding',
                  'lncRNA': 'lncRNA',
                  'pseudogene': 'Pseudogene'}
        # biotypes.reverse()
        c_dict, order = get_shade_colors(c, biotypes)
        # order.reverse()

    df = df.pivot(index=feat_col, columns=['biotype_category'])
    df.columns = df.columns.droplevel(0)
    df.reset_index(inplace=True)

    df[feat_col] = df[feat_col].astype(str)
    x = df[feat_col].unique().tolist()

    # loop through biotypes
    bottom = [0 for i in range(len(x))]
    for b in order:
        y = df[b].tolist()
        plt.bar(x, y, color=c_dict[b], bottom=bottom)
        bottom = [b_coord+y_coord for b_coord, y_coord in zip(bottom, y)]

    plt.xlabel('# transcripts / gene')
    plt.ylabel('# genes')
    sns.despine()


    leg_labels = [b_dict[o] for o in order]
    plt.legend(leg_labels, bbox_to_anchor=(.6, 1.05))
    ax = plt.gca()
    leg = ax.get_legend()
    shade_dict, _ = get_shade_colors('#000000', order)
    for i, o in enumerate(order):
        leg.legendHandles[i].set_color(shade_dict[o])

    # ax.tick_params(axis="x", rotation=45)


    fname = f'{opref}/isos_per_gene_support.png'
    print(fname)
    plt.savefig(fname, dpi=500, bbox_inches='tight')

    fname = f'{opref}/isos_per_gene_support.pdf'
    plt.savefig(fname, dpi=500, bbox_inches='tight')

    return df

def add_bool_heatmap(tpm_df, ax, species):
    if species == 'human':
        keep_cols = ['triplet', 'gene_name', 'is_mane_orf',
                     'full_orf', 'nmd', 'Known']
    else:
         keep_cols = ['triplet', 'gene_name',
                     'full_orf', 'nmd', 'Known']
    temp = tpm_df[keep_cols]
    temp['iso_trip'] = temp.gene_name+' '+temp.triplet
    temp.drop(['triplet', 'gene_name'], axis=1, inplace=True)
    temp.set_index('iso_trip', inplace=True)
    temp.index.name = ''
    if species == 'human':
        r_map = {'is_mane_orf': 'MANE ORF',
                 'full_orf': 'Full ORF',
                 'nmd': 'NMD'}
    else:
        r_map = {'full_orf': 'Full ORF',
                 'nmd': 'NMD'}
    temp.rename(r_map, axis=1, inplace=True)
    if species == 'human':
        bool_cols = ['Known', 'MANE ORF', 'Full ORF', 'NMD']
    else:
        bool_cols = ['Known', 'Full ORF', 'NMD']

    temp = temp[bool_cols]
    temp2 = pd.DataFrame()

    m = {True: '*', False: ''}
    for c in bool_cols:
        temp2[c] = temp[c].map(m)


    # if species == 'mouse':
    #     temp.drop('MANE ORF', axis=1, inplace=True)
    #     temp2.drop('MANE ORF', axis=1, inplace=True)
    ax = sns.heatmap(temp, cbar=False, cmap='Purples',
                     linewidths=0.5, linecolor='k',
                     annot=temp2, square=True,
                     fmt='', ax=ax)
    ax.tick_params(left=False,
                   right=False, labelright=False,
                   bottom=False, labelbottom=False,
                   labeltop=True)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    return ax

def get_cds(exons, orf_start, orf_end, strand):
    # ChatGPT wrote most of this thx friend

    def rev_exons(exons):
        return [(e[1],e[0]) for e in exons[::-1]]

    if strand == '-':
        exons = rev_exons(exons)
        orf_start, orf_end = orf_end, orf_start

    cds = []
    for exon_start, exon_end in exons:

        # determine the overlap between the current exon and the ORF
        overlap_start = max(exon_start, orf_start)
        overlap_end = min(exon_end, orf_end)

        # if there is no overlap, skip this exon
        if overlap_start > overlap_end:
            continue

        # add the CDS range for this exon to the list
        cds.append((overlap_start, overlap_end))

    if strand == '-':
        cds = rev_exons(cds)

    return cds

def plot_browser_isos(ca, sg, gene,
                      obs_col, obs_condition,
                      filt_ab, pp_summary,
                      major_set,
                      h=0.1, w=56, fig_w=14, species='human',
                      add_tss=False,
                      add_ccre=False,
                      major=False,
                      order='expression',
                      light_shade=None,
                      dark_shade=None):
    """
    Plot browser style isoform models for a given sample
    """

    def plot_tss(ca, sg, tpm_df, x, y, h, ax):
        tpm_df = add_feat(tpm_df, kind='tss', col='transcript_id')
        tpm_df.head()
        tss_df = ca.tss.loc[ca.tss.Name.isin(tpm_df.tss.tolist())]
        tss_df.head()
        regions = [(entry.Start, entry.End) for ind, entry in tss_df.iterrows()]
        color = get_sector_colors()[0]['tss']
        ax = sg.pg.plot_regions(regions, sg.pg.scale, sg.pg.strand, sg.pg.g_min, sg.pg.g_max, x, y, h, color, ax)
        return ax

    def plot_cds(entry, sg, x, y, h, ax, color):
        c_dict, order = get_sector_colors()

        def get_inc_region(entry, how, scale):
            inc = 80
            if how == 'start' and entry.Strand == '+':
                col = 'CDS_Start'
            elif how == 'stop' and entry.Strand == '+':
                col = 'CDS_Stop'
            elif how == 'start' and entry.Strand == '-':
                col = 'CDS_Stop'
            elif how == 'stop' and entry.Strand == '-':
                col = 'CDS_Start'
            if entry.Strand == '+':
                regions = (entry[col],entry[col]+inc)
            elif entry.Strand == '-':
                regions = (entry[col],entry[col]-inc)
            regions = [regions]
            # print(f'Scale: {scale}')
            # print(f'CDS inc: {inc}')

            return regions

        # get start and stop
        strand = entry.Strand
        if strand == '+':
            orf_start = entry['CDS_Start']
            orf_end = entry['CDS_Stop']
        elif strand == '-':
            orf_end = entry['CDS_Start']
            orf_start = entry['CDS_Stop']
        loc_path = sg.pg.loc_path
        exons = [(sg.pg.loc_df.loc[v1, 'coord'],
            sg.pg.loc_df.loc[v2, 'coord']) \
            for v1,v2 in zip(loc_path[:-1],loc_path[1:])][::2]
        cds = get_cds(exons, orf_start, orf_end, strand)
        ax = sg.pg.plot_regions(cds, sg.pg.scale, sg.pg.strand, sg.pg.g_min, sg.pg.g_max, x, y, h, color, ax)

        # start codon
        color = c_dict['tss']
        regions = get_inc_region(entry, 'start', sg.pg.scale)
        ax = sg.pg.plot_regions(regions, sg.pg.scale, sg.pg.strand, sg.pg.g_min, sg.pg.g_max, x, y, h, color, ax)

        # stop codon
        color = c_dict['tes']
        regions = get_inc_region(entry, 'stop', sg.pg.scale)
        ax = sg.pg.plot_regions(regions, sg.pg.scale, sg.pg.strand, sg.pg.g_min, sg.pg.g_max, x, y, h, color, ax)

        return ax


    def plot_ccre(ca, sg, x, y, h, ax):

        # ccre regions
        sources = ['pls', 'pels', 'dels']
        ccre = ca.tss_map.loc[ca.tss_map.source.isin(sources)].copy(deep=True)

        # get ranges w/i this region
        min_coord = sg.pg.g_min
        max_coord = sg.pg.g_max
        chrom = sg.pg.chrom

        ccre['min_coord'] = ccre[['Start', 'End']].min(axis=1)
        ccre['max_coord'] = ccre[['Start', 'End']].max(axis=1)

        # subset on regions
        ccre = pr.PyRanges(ccre)
        region = pr.from_dict({'Chromosome': [chrom],
                                        'Start': [min_coord],
                                        'End': [max_coord]})
        ccre = ccre.intersect(region, strandedness=None)
        ccre = ccre.as_df()

        # colors
        c_dict, _ = get_ccre_colors()
        colors = [c_dict[s] for s in ccre.source.tolist()]
        regions = [(entry.Start, entry.End) for ind, entry in ccre.iterrows()]
        ax = sg.pg.plot_regions(regions, sg.pg.scale, sg.pg.strand, sg.pg.g_min, sg.pg.g_max, x, y, h, colors, ax)

        return ax

    def get_major_isos(major_set, gene, sample=None):
        """
        Get list of major isfoorms in a given sample
        """
        df = pd.read_csv(major_set, sep='\t')
        df = df.loc[df.gname == gene]
        if sample:
            df = df.loc[df['sample'] == sample]
        tids = df.tid.unique().tolist()
        return tids

    def get_isos(ca, filt_ab, gene, sample, species):
        df = pd.read_csv(filt_ab, sep='\t')
        df = get_det_table(df,
                       groupby='sample',
                       how='iso',
                       min_tpm=1,
                       gene_subset='polya',
                       species=species)
        df = df.loc[sample]
        df = df.to_frame()
        df = df.loc[df[sample]==True]
        gid = ca.triplets.loc[ca.triplets.gname==gene, 'gid'].values[0]
        df.reset_index(inplace=True)
        df['gid'] = df['index'].str.split('[', expand=True)[0]
        df = df.loc[df.gid == gid]
        tids = df['index'].tolist()
        return tids


    def get_tpm_df(sg, tids, obs_col, obs_condition):
        tpm_df = swan.calc_tpm(sg.adata, obs_col=obs_col).sparse.to_dense()
        tpm_df = tpm_df.transpose()
        tpm_df = tpm_df.loc[tids, obs_condition].to_frame()
        return tpm_df

    if major:
        tids = get_major_isos(major_set, gene, obs_condition)
    else:
        tids = get_isos(ca, filt_ab, gene, obs_condition, species)
    tpm_df = get_tpm_df(sg, tids, obs_col, obs_condition)

    # colormap definition
    if not light_shade:
        light_shade = get_sector_colors()[0]['mixed']
    if not dark_shade:
        dark_shade = get_sector_colors()[0]['simple']
    cmap = mpl.colors.LinearSegmentedColormap.from_list('', [light_shade, dark_shade])

    # figure / subplot settings
    fig_h = len(tids)/4
    w1_rat = (fig_w/6)*2
    w2_rat = (fig_w/6)*4
    fig, (ax, ax2) = plt.subplots(1,2, figsize=(fig_w, fig_h),
                                  gridspec_kw={'width_ratios':(w1_rat, w2_rat)},
                                  frameon=False, sharey=True)
    fig.subplots_adjust(wspace=0.00)
    fig.subplots_adjust(hspace=0.00)

    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    # plotting order
    tpm_df = tpm_df.sort_values(by=obs_condition, ascending=False)
    tpm_df = tpm_df[[obs_condition]]

    if order == 'tss':
        tpm_df = add_feat(tpm_df, kind='tss', col='index')
        tpm_df.sort_values(by=['tss', obs_condition], ascending=[True, False], inplace=True)
        tpm_df.drop('tss', axis=1, inplace=True)

    # add reference ids
    if species == 'human':
        ref_source = 'v40'
    elif species == 'mouse':
        ref_source = 'vM25'
    df = ca.t_map.loc[ca.t_map.source == ref_source]
    df = df[['transcript_id','original_transcript_id', 'original_transcript_name']]
    tpm_df = tpm_df.merge(df, how='left', left_index=True, right_on='transcript_id')

    # clean up for transcripts that aren't known
    df = ca.t_map.loc[ca.t_map.source == 'lapa']
    df = df[['gene_name', 'transcript_id', 'transcript_name',
             'tss_id', 'ic_id', 'tes_id']].drop_duplicates()
    tpm_df = tpm_df.merge(df, how='left', on='transcript_id')
    tpm_df.fillna('N/A', inplace=True)

    if species == 'human':
        refs = ['v40', 'v29']
    elif species == 'mouse':
        refs = ['vM21', 'vM25']
    known_tids = ca.t_map.loc[ca.t_map.source.isin(refs)].transcript_id.unique().tolist()
    tpm_df['Known'] = tpm_df.transcript_id.isin(known_tids)

    # triplets rather than entire transcript name
    tpm_df['triplet'] = tpm_df.transcript_id.str.split('[', n=1, expand=True)[1]
    tpm_df['triplet'] = tpm_df.triplet.str.split(']', n=1, expand=True)[0]
    tpm_df['triplet'] = '['+tpm_df.triplet+']'

    # NMD and ORF
    # TODO - update
    print('PLEASE UPDATE ME WHEN U HAVE PROTEIN RESULTS')
    if pp_summary:
        pp_df = pd.read_csv(pp_summary, sep='\t')
        pp_df.rename({'tid': 'transcript_id'}, axis=1, inplace=True)
        tpm_df = tpm_df.merge(pp_df[['transcript_id', 'nmd', 'full_orf',
                                     'seq', 'Strand', 'CDS_Start', 'CDS_Stop']],
                              how='left',
                              on='transcript_id')

        # are these the mane orf
        if species == 'human':
            tpm_df['gid'] = tpm_df['ic_id'].str.split('_', expand=True)[0]
            gid = tpm_df.gid.values[0]
            mane_orf = get_mane_orf(pp_summary,
                                    'v40_cerberus',
                                    gid=gid).seq.values[0]
            tpm_df['is_mane_orf'] = tpm_df['seq'] == mane_orf

    else:
        tpm_df[['nmd', 'full_orf', 'is_mane_orf',
                'seq', 'Strand', 'CDS_Start', 'CDS_Stop']] = False
        tpm_df[['CDS_Start', 'CDS_Stop']] = 0

    # first add the labels
    ax = add_bool_heatmap(tpm_df, ax, species)

    # get y coords from heatmap
    y_locs = ax.get_yticks()

    # get x coords from heatmap
    x = 0
    i = 0

    # get height of cdss
    h_cds = h*2

    linewidth = 1.5

    for index, entry in tpm_df.iterrows():

        y_ytick = y_locs[i]

        # y coords
        y = y_ytick-(h/2)
        y_cds = y_ytick-(h_cds/2)

        # tid
        tid = entry['transcript_id']

        # color by TPM
        if len(tpm_df.index) == 1:
            norm_val = entry[obs_condition]
        else:
            norm_val = (entry[obs_condition]-tpm_df[obs_condition].min())/(tpm_df[obs_condition].max()-tpm_df[obs_condition].min())
        color = cmap(norm_val)

        # plot models + cds
        ax2 = sg.plot_browser(tid, y=y, x=x, h=h, w=w, color=color,
                              ax=ax2, linewidth=linewidth)

        print('TODO - also add this back in')
        # ax2 = plot_cds(entry, sg, x, y_cds, h_cds, ax2, color)

        i+=1

    y_space=0.5

    i = 1
    if add_tss:
        y_ytick = (len(tpm_df.index) + i)
        y = y_ytick-(h_cds/2)
        ax2 = plot_tss(ca, sg, tpm_df, x, y, h_cds, ax2)
        i += 1

    if add_ccre:
        y_ytick = (len(tpm_df.index) + i)
        y = y_ytick-(h_cds/2)
        ax2 = plot_ccre(ca, sg, x, y, h_cds, ax2)
        i += 1

    # add scale
    y_ytick = (len(tpm_df.index) + i)
    y = y_ytick-(h_cds/2)
    ax2 = sg.pg.plot_scale(x, y, h_cds, w, ax2, linewidth=linewidth)

    # print(sg.pg.scale*(250000))

    plt.tight_layout()
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    ax2.set_frame_on(False)

    return ax, tpm_df

def plot_browser_isos_2(h5,
                        swan_file,
                        filt_ab,
                        pp_summary,
                        major_isos,
                        gene,
                        obs_col,
                        obs_condition,
                        ofile,
                        **kwargs):

    ca = cerberus.read(h5)
    sg = swan.read(swan_file)

    ax, tpm_df = plot_browser_isos(ca,
                                   sg,
                                   gene,
                                   obs_col,
                                   obs_condition,
                                   filt_ab,
                                   pp_summary,
                                   major_isos,
                                   h=0.4,
                                   w=10,
                                   fig_w=6,
                                   **kwargs)
    plt.savefig(ofile, dpi=500, bbox_inches='tight')

    return ax, tpm_df

def plot_gene_tpm_v_predom_t_pi(h5,
                                major_isos,
                                source,
                                obs_col,
                                obs_condition,
                                gene_subset,
                                ofile,
                                label_genes=None):
    ca = cerberus.read(h5)
    df = ca.triplets.copy(deep=True)
    df = df.loc[df.source == source]
    df = df.loc[df[obs_col] == obs_condition]
    sns.set_context('paper', font_scale=2)


    # get predominant transcripts in this sample
    prin_isos = pd.read_csv(major_isos, sep='\t')
    prin_isos = prin_isos.loc[prin_isos[obs_col] == obs_condition]
    prin_isos = prin_isos.loc[prin_isos.pi_rank == 1]
    prin_isos['gid_stable'] = cerberus.get_stable_gid(prin_isos, 'gid')
    prin_isos.drop(['gid', 'gname'], axis=1, inplace=True)
    prin_isos.rename({'gid_stable':'gid'}, axis=1, inplace=True)

    import pdb; pdb.set_trace()

    df = df.merge(prin_isos, how='left', on=['gid', 'sample'])
    df['log2tpm'] = np.log2(df.gene_tpm+1)

    # genes with one isoform
    df['one_iso'] = df.n_iso==1

    # sample and one iso color
    c_dict, order = get_biosample_colors()
    color = c_dict[obs_condition]
    c_dict = {True: '#b7b7b7', False: color}

    x_col = 'pi'
    y_col = 'log2tpm'

    pi_dict = {'low': 50, 'high': 90}
    tpm_dict = {'low': 20, 'high': 100}

    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    ax = sns.jointplot(data=df, x=x_col, y=y_col,
                       hue='one_iso',
                       palette=c_dict, size=1, alpha=0.5,
                       marginal_kws={'linewidth':0})
    ax = ax.ax_joint
    xlabel = f'Pi of predominant transcript in {obs_condition}'
    ylabel = f'log2(TPM+1) of gene in {obs_condition}'
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.get_legend().remove()

    if label_genes:
        for g in label_genes:
            x = df.loc[df.gname == g, x_col].values[0]
            y = df.loc[df.gname == g, y_col].values[0]
            ax.annotate(g, (x,y), fontsize='small', fontstyle='italic',
                        xytext=(4,-4.5), textcoords='offset pixels',
                        arrowprops={'width':10, 'headwidth':0})

        # add sector lines
        xlims = ax.get_xlim()
        ylims = ax.get_ylim()

        # tpm
        low_y = math.log2(tpm_dict['low']+1)
        high_y = math.log2(tpm_dict['high']+1)

        # pi
        low_x = pi_dict['low']
        high_x = pi_dict['high']
        color = '#5c5c5c'
        ax.hlines(low_y, xlims[0], xlims[1],
                      colors=color, linestyles='dashed',
                      linewidth=2)
        ax.hlines(high_y, xlims[0], xlims[1],
                      colors=color, linestyles='dashed',
                      linewidth=2)

        ax.vlines(low_x, ylims[0], ylims[1],
                      colors=color, linestyles='dashed',
                      linewidth=2)
        ax.vlines(high_x, ylims[0], ylims[1],
                      colors=color, linestyles='dashed',
                      linewidth=2)

        plt.savefig(ofile, dpi=500)

        # quantify how many genes fall within certain thresholds
        pi_col = 'pi'
        tpm_col = 'gene_tpm'
        n = len(df.index)
        for pi_key, pi in pi_dict.items():
            for tpm_key, tpm in tpm_dict.items():
                print('{} pi and {} tpm'.format(pi_key, tpm_key))
                if pi_key == 'low':
                    temp = df.loc[df[pi_col]<pi].copy(deep=True)
                    pi_label = '<'
                elif pi_key == 'high':
                    temp = df.loc[df[pi_col]>pi].copy(deep=True)
                    pi_label = '>'
                if tpm_key == 'low':
                    temp = temp.loc[temp[tpm_col]<tpm]
                    tpm_label = '<'
                elif tpm_key == 'high':
                    temp = temp.loc[temp[tpm_col]>tpm]
                    tpm_label = '>'
                n_num = len(temp.index)
                print('{:.2f}% ({}/{}) of protein coding genes in ovary have pi {} {} and tpm {} {}'.format((n_num/n)*100, n_num, n, pi_label, pi, tpm_label, tpm))
                print()

        # what about how many are w/i pi thresholds but over tpm threshold?
    #     pi_dict = {'low': 50, 'high': 90}
    # tpm_dict = {'low': 20, 'high': 100}
        tpm_min = tpm_dict['high']
        temp = df.loc[df[tpm_col]>tpm_min]
        n = len(temp.index)
        for pi_key, pi in pi_dict.items():
            if pi_key == 'low':
                    temp2 = temp.loc[temp[pi_col]<pi].copy(deep=True)
                    pi_label = '<'
            elif pi_key == 'high':
                temp2 = temp.loc[temp[pi_col]>pi].copy(deep=True)
                pi_label = '>'
            n_num = len(temp2.index)
            print('{:.2f}% ({}/{}) of protein coding genes >{} tpm in ovary have pi {} {}'.format((n_num/n)*100, n_num, n, tpm_min, pi_label, pi))
            print()

def plot_n_predom_transcripts(pi_tpm_file,
                              filt_ab,
                              ver,
                              gene_subset,
                              min_tpm,
                              fname,
                              obs_col='sample',
                              species='human',
                              max_isos=None,
                              figsize=(6,6)):
    df = pd.read_csv(pi_tpm_file, sep='\t')

    ab_df = pd.read_csv(filt_ab, sep='\t')
    det_df = get_det_table(ab_df,
               how='iso',
               min_tpm=min_tpm,
               gene_subset=gene_subset,
               groupby=obs_col,
               species=species)

    # only predominant transcripts
    df = df.loc[df.triplet_rank==1]

    # only expressed transcripts in each sample
    det_df = det_df.melt(ignore_index=False,
                var_name='tid',
                value_name='det').reset_index()
    if obs_col=='sample':
        det_df.rename({'biosample': 'sample'}, axis=1)

    det_df = det_df.loc[det_df.det==True]
    print(len(df.index))
    df = df.merge(det_df, how='inner', on=[obs_col, 'tid'])
    print(len(df.index))


    # count number of unique predominant transcripts
    df_back = df.copy(deep=True)
    df = df[['tid', 'gid']].groupby(['gid']).nunique().reset_index()
    df.rename({'tid': 'n_predom_ts'}, axis=1, inplace=True)

    # limit to gene subset
    if gene_subset:
        gene_df, _, _ = get_gtf_info(how='gene',
                                     ver=ver,
                                     add_stable_gid=True)
        gene_df = gene_df[['gid_stable', 'biotype', 'gname']]
        df = df.merge(gene_df, how='left',
                      left_on='gid', right_on='gid_stable')
        df = df.loc[df.biotype==gene_subset]
        df.drop(['biotype', 'gid_stable'], axis=1, inplace=True)

    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    color = get_talon_nov_colors(['Known'])[0]['Known']

    # a wee bit of math
    print(f'Median predominant transcripts / gene: {df.n_predom_ts.median()}')
    df['one_iso'] = df.n_predom_ts == 1
    temp = df[['n_predom_ts', 'one_iso']].groupby('one_iso').count().reset_index()
    n = temp.n_predom_ts.sum(axis=0)
    n_num = temp.loc[temp.one_iso == False, 'n_predom_ts'].values[0]
    print(n)
    print(n_num)
    print(f'{n_num}/{n} {(n_num/n)*100:.2f}% protein-coding genes have >1 predominant isoforms across samples')


    if max_isos:
        df.loc[df.n_predom_ts>=max_isos, 'n_predom_ts'] = max_isos
        xticks = [i for i in range(0, max_isos+1, 5)]
        xtick_labels = [str(xtick) for xtick in xticks]
        xtick_labels[-1] = f'{max_isos}+'


    # and make a beautiful plot
    plt.figure(figsize=figsize)
    height = figsize[1]
    width = figsize[0]
    aspect = width/height

    ax = sns.displot(df, x='n_predom_ts', kind='hist',
             discrete=True,
             color=color,
             linewidth=0,
             alpha=1,
             height=height,
                     aspect=aspect)
    xlabel = '# predominant transcripts'
    ylabel = '# genes'
    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    if max_isos:
        sub_ax = plt.gca()
        sub_ax.set_xticks(xticks)
        sub_ax.set_xticklabels(xtick_labels)

    plt.savefig(fname, dpi=500, bbox_inches='tight')

    # ax.tick_params(axis="x", rotation=90)

    return df_back, df

def make_triplet_feat_upset(h5,
                            filt_ab,
                            feat,
                            gene_subset,
                            min_tpm,
                            ofile):

    ca = cerberus.read(h5)
    ids = get_det_feats(h5,
                    filt_ab,
                    feat,
                    how=feat,
                    gene_subset=gene_subset,
                    min_tpm=min_tpm)

    if feat == 'tss':
        df = ca.tss.copy(deep=True)
        m = {'GENCODE': ['v29', 'v40'],
             'CAGE, RAMPAGE': ['encode_cage', 'encode_rampage', 'fantom_cage'],
             'GTEx': ['gtex'],
             'PLS cCREs': ['pls'],
             'pELS/dELS cCREs': ['pels', 'dels']}
             # 'cCREs': ['pls', 'dels', 'pels']}
    elif feat == 'ic':
        df = ca.ic.copy(deep=True)
        m = {'GENCODE': ['v29', 'v40'],
             'GTEx': ['gtex']}
    elif feat == 'tes':
        df = ca.tes.copy(deep=True)
        m = {'GENCODE': ['v29', 'v40'],
             'PAS-seq, PolyA Atlas': ['pas', 'polya_atlas'],
             'GTEx': ['gtex']}

    print(feat)
    print(m)

    # subset on detected feats
    df = df.loc[df.Name.isin(ids)]

    # get melted version
    end_upset = upsetplot.from_memberships(df.source.str.split(','), data=df)
    all_sources = end_upset.index.names
    end_upset.reset_index(inplace=True)
    sources = []
    for key, item in m.items():
             end_upset[key] = end_upset[item].any(axis=1)
             sources.append(key)
    end_upset.drop(all_sources, axis=1, inplace=True)
    end_upset.set_index(sources, inplace=True)

    # # filter for given sources
    #
    # sources = list(set(all_sources) - set(['lapa']))
    # end_upset.reset_index(inplace=True)
    # end_upset.drop('lapa', axis=1, inplace=True)
    # end_upset.set_index(sources, inplace=True)

    # # limit to number of subsets
    # if n_subset:
    #     uniques, counts = np.unique(end_upset.index, return_counts=True)
    #     sorted_uniques = [x for _, x in sorted(zip(counts, uniques), reverse=True)]
    #     end_upset = end_upset.loc[sorted_uniques[:n_subset]].copy(deep=True)

    # make the plot
    c_dict, _ = get_feat_colors()
    c = c_dict[feat]
    # tss should be longer
    if feat == 'tss':
        fig = plt.figure(figsize=(22,5))
    else:
        fig = plt.figure(figsize=(11,5))
    sns.set_context('paper', font_scale=1.8)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    upsetplot.plot(end_upset, subset_size='auto',
                    show_counts='%d', sort_by='cardinality',
                    facecolor=c, fig=fig, shading_color='white', element_size=None)
    plt.savefig(ofile, dpi=500, bbox_inches='tight')


def plot_feat_len_hist(cerberus_h5,
                       filt_ab,
                       feat,
                       gene_subset,
                       min_tpm,
                       ofile):
    ca = cerberus.read(cerberus_h5)

    ids = get_det_feats(cerberus_h5,
                        filt_ab,
                        feat,
                        how=feat,
                        gene_subset=gene_subset,
                        min_tpm=min_tpm)
    if feat == 'tss':
        df = ca.tss.loc[ca.tss.Name.isin(ids)]
    elif feat == 'tes':
        df = ca.tes.loc[ca.tes.Name.isin(ids)]
    df['region_len'] = abs(df.Start-df.End)
    print(df.region_len.min())

    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    c_dict, order = get_feat_colors(feat)
    color = c_dict[feat]

    ax = sns.displot(df, x='region_len', kind='hist',
             linewidth=0,
             color=color,
             alpha=1,
             binwidth=25,
             edgecolor=None,
             log_scale=(False, True))


    ax = plt.gca()
    ylim = ax.get_ylim()
    ax.set_ylim(10**-.3, ylim[1])

    df.region_len.max()
    ylabel = f'# detected {feat.upper()}s'
    xlabel = f'{feat.upper()} length (bp)'

    ax.set(xlabel=xlabel, ylabel=ylabel)
    plt.savefig(ofile, dpi=700, bbox_inches='tight')

    # some math
    bounds_list = [(0,500), (250,df.region_len.max())]
    for bounds in bounds_list:
        n_num = len(df.loc[(df.region_len>=bounds[0])&(df.region_len<=bounds[1])].index)
        n = len(df.index)
        print(f'{(n_num/n)*100:.2f}% of regions ({n_num}/{n}) are b/w {bounds[0]} and {bounds[1]} bp long')

    return df


def plot_ic_novelty(fname,
                    source,
                    oprefix,
                    ylim=None,
                    pass_list=None,
                    novs=None,
                    save_type='pdf'):
    """
    Plot number of intron chains per novelty category.

    Parameters:
        fname (str): Cerberus annotation file name
        source (str): Source in cerberus annotation
        oprefix (str): Place to save
        ylim (int): y limit of resultant plot
        pass_list (list of str): List of ic IDs to retain
        novs (list of str): Novelty types to include
        save_type (str): Choose from 'pdf' or 'png'
    """

    # sns.set_context('paper', font_scale=1.6)

    ca = cerberus.read(fname)


    temp = ca.t_map.loc[ca.t_map.source==source].copy(deep=True)
    temp = temp[['ic_id']]
    nov = ca.ic[['Name', 'novelty']]
    temp = temp.merge(nov, how='left', left_on='ic_id', right_on='Name')
    temp.drop_duplicates(inplace=True)

    if pass_list:
        temp = temp.loc[temp.ic_id.isin(pass_list)]

    temp = temp[['ic_id', 'novelty']]
    temp = temp.groupby('novelty').count()

    temp.reset_index(inplace=True)
    temp.rename({'ic_id': 'counts'}, axis=1, inplace=True)
    print(temp)
    if not novs:
        novs = temp.novelty.unique().tolist()
    c_dict, order = get_ic_nov_colors(cats=novs)

    temp = temp.loc[temp.novelty.isin(novs)].copy(deep=True)
    complete = temp[['counts']].sum(axis=0)
    print('Number of complete intron chains: {}'.format(complete))

    # actual plotting
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    # plt.figure(figsize=(4,6))
    plt.figure(figsize=(3,4))

    g = sns.catplot(data=temp, x='novelty',
                y='counts', kind='bar',
                saturation=1,
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=45) for ax in g.axes.flat]
    g.set_ylabels('# intron chains')
    g.set_xlabels('')

    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, temp, 'counts')

    if ylim:
        g.set(ylim=(0,ylim))

    # save figure
    fname = '{}_ic_novelty'.format(oprefix)
    g.savefig(fname+'.png', dpi=500, bbox_inches='tight')
    g.savefig(fname+'.pdf', dpi=500, bbox_inches='tight')

    plt.show()
    plt.clf()

    def plot_exp_v_iso_biotype_boxplot(h5,
                                       ver,
                                       ofile):
        ca = cerberus.read(h5)

        # limit to the sample dets
        df = ca.triplets.loc[ca.triplets.source=='sample_det'].copy(deep=True)

        # tpm bins
        tpm_bins = [1, 10, 100, df.gene_tpm.max()]
        labels = ['Low (1-10)', 'Medium (10-100)', 'High (100-max)']

        # get the most highly-expressed sample per gene
        df = df.sort_values(by=['gid', 'gene_tpm'], ascending=[False,False])
        df = df.drop_duplicates(subset='gid', keep='first')

        # group into bins
        df['tpm_bin'] = pd.cut(df.gene_tpm, tpm_bins, labels=labels)
        df['other_tpm_bin'] = pd.cut(df.gene_tpm, tpm_bins)
        print(df.other_tpm_bin.unique())

        # add gene biotypes
        gene_df, _, _ = get_gtf_info(how='gene', ver=ver)
        gene_df['gid'] = cerberus.get_stable_gid(gene_df, col='gid')
        df = df.merge(gene_df[['gid', 'biotype_category']], how='left', on='gid')

        order = get_polya_cats()
        disp_dict = {'protein_coding': 'Protein coding',
                     'lncRNA': 'lncRNA',
                     'pseudogene': 'Pseudogene'}
        order = [disp_dict[b] for b in order]
        df['biotype_category_disp'] = df.biotype_category.map(disp_dict)

        c = get_talon_nov_colors()[0]['Known']
        c_dict, order = get_shade_colors(c, order)

        sns.set_context('paper', font_scale=1.8)
        mpl.rcParams['font.family'] = 'Arial'
        mpl.rcParams['pdf.fonttype'] = 42

        ax = sns.catplot(data=df,
                         x='tpm_bin',
                         y='n_iso',
                         hue='biotype_category_disp',
                         kind='box',
                         hue_order=order,
                         palette = c_dict,
                         fliersize=1,
                         linewidth=1,
                         height=4, aspect=(5/4),
                         saturation=1)

        ax.set(ylim=(0,50), xlabel='Max. gene expression (TPM)', ylabel='# transcripts / gene')
        ax.tick_params(axis="x", rotation=45)

        fname = ofile
        plt.savefig(fname, dpi=500, bbox_inches='tight')

        return df

def plot_ends_per_ic(df, ca,
                     feat,
                     fname,
                     rm_monoexonic=True):
    """
    Plot a histogram of the # ends (tss or tes) per ic

    Parameters:
        df (pandas DataFrame): DF where index is transcript
            id w/ triplet form
        feat (str): {'tss', 'tes'}
        fname (str): Output file name to save to

    Returns:
        temp (pandas DataFrame): DF w/ # ends / ic
    """

    # plot settings
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    sns.set_context('paper', font_scale=1.8)


    temp = df.copy(deep=True)
    feats = ['ic', feat]
    for feat in feats:
        temp = add_feat(temp, col='index', kind=feat)

    # keep only the relevant features and drop duplicated
    # combinations
    temp.reset_index(drop=True, inplace=True)
    temp = temp[[feat, 'ic']].copy(deep=True)
    temp.drop_duplicates(inplace=True)

    # groupby and count number of feature per ic
    temp = temp.groupby('ic').nunique().reset_index()

    # merge with novelty info from h5
    temp = temp.merge(ca.ic[['Name', 'Coordinates']],
                      how='left', left_on='ic', right_on='Name')

    if rm_monoexonic:
        temp = temp.loc[temp.Coordinates != '-']
    temp.drop(['Name', 'Coordinates'], axis=1, inplace=True)

    # plotting
    sns.set_context('paper', font_scale=1.8)
    c_dict, order = get_end_colors()
    c = c_dict[feat]
    # ax = sns.displot(temp, x=feat, kind='hist',
    #                  linewidth=0,
    #                  color=c,
    #                  discrete=True,
    #                  alpha=1)
    ax = sns.displot(temp, x=feat, kind='hist',
                 linewidth=0,
                 color=c,
                 discrete=True,
                 alpha=1,
                 binwidth=1,
                 log_scale=(False, True))


    ylabel = '# ICs'
    xlabel = '# {}s / IC'.format(feat.upper())
    # ax.fig.get_axes()[0].set_yscale('log')
    # ax.fig.get_axes()[0].xaxis.set_major_formatter(FormatStrFormatter('%d'))
    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    ax = plt.gca()
    ylim = ax.get_ylim()
    ax.set_ylim(10**-.3, ylim[1])

    # _ = ax.set_yscale("log")

    # plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.savefig(fname, dpi=700)

    # do a lil math
    temp['one_end'] = temp[feat]==1
    temp2 = temp[['ic', 'one_end']].groupby('one_end').count().reset_index()
    n = temp2.ic.sum(axis=0)
    n_num = temp2.loc[temp2.one_end == True, 'ic'].values[0]
    print(f'{(n_num/n)*100:.2f}% ({n_num}/{n} of unique ics have 1 {feat}')

    return temp
