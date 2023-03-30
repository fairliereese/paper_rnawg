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


from .utils import *

def get_ccre_colors():
    pls = '#FF0000'
    pels = '#FFA700'
    dels = '#FFCD00'
    order = ['pls', 'pels', 'dels']
    c_dict = {'pls': pls,
              'pels': pels,
              'dels': dels}
    return c_dict, order

def get_not_det_color():
    return '#E5ECF6'

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
    fname = '{}/../mouse/refs/tissue_colors.tsv'.format(d)
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
    c_dict['adrenal'] = '#8e361b'
    c_dict['hippocampus'] = '#9a3c4f'
    c_dict['cortex'] = '#634273'
    
    order += ['adrenal', 'hippocampus', 'cortex']
    

    # f1219
    c_dict['f1219'] = '#4340bc'
    order += ['f1219']

    return c_dict, order

def get_shade_colors(color, order):
    c_dict = {}
    min_color = '#FFFFFF'
    cmap = mpl.colors.LinearSegmentedColormap.from_list('temp', [color, min_color], N=len(order)+1)
    for i, cat in enumerate(order):
        c_dict[cat] = mpl.colors.to_hex(cmap(i))

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

def get_tissue_cell_line_colors(cats=None):
    tissue = '#e39f24'
    cell_line = '#7680e8'
    c_dict = {'cell_line': cell_line,
              'tissue': tissue}
    order = ['cell_line', 'tissue']
    c_dict, order = rm_color_cats(c_dict, order, cats)
    return c_dict, order

def get_support_sector_colors(sector=None):
    c_dict, order = get_sector_colors()

    # get different color for each age / tissue
    new_c_dict = {}
    min_color = '#FFFFFF'
    ages = [False, True]
    order = []
    for tissue, max_color in c_dict.items():
        cmap = mpl.colors.LinearSegmentedColormap.from_list(tissue, [min_color, max_color], N=len(ages)+1)
        for i, age in enumerate(ages):
            key = '{}_{}'.format(tissue, age)
            new_c_dict[key] = mpl.colors.to_hex(cmap(i+1))
            order.append(key)

    if sector:
        cats = [key for key in new_c_dict.keys() if sector in key]
        new_c_dict, order = rm_color_cats(new_c_dict, order, cats)
        new_c_dict['{}_False'.format(sector)]

        c = dict()
        o = []
        for b in [True, False]:
            o.append(b)
            key = '{}_{}'.format(sector, b)
            c[b] = new_c_dict[key]
        new_c_dict = c
        order = o

    return new_c_dict, order

def get_feat_triplet_colors(cats=None):
    tss = '#56B4E9'
    tes = '#E69F00'
    splicing = '#CC79A7'
    triplet = '#009E73'
    c_dict = {'tss': tss,
              'ic': splicing,
              'tes': tes,
              'triplet': triplet}
    order = ['triplet', 'tss', 'ic', 'tes']
    
    c_dict, order = rm_color_cats(c_dict, order, cats)
    return c_dict, order

def get_feat_triplet_colors_2(cats=None):
    tss = '#56B4E9'
    tes = '#E69F00'
    splicing = '#CC79A7'
    triplet = '#009E73'
    c_dict = {'tss': tss,
              'ic': splicing,
              'tes': tes,
              'iso': triplet}
    order = ['iso', 'tss', 'ic', 'tes']
    
    c_dict, order = rm_color_cats(c_dict, order, cats)
    return c_dict, order

def get_feat_colors(cats=None):
    tss = '#56B4E9'
    tes = '#E69F00'
    splicing = '#CC79A7'
    c_dict = {'tss': tss,
              'ic': splicing,
              'tes': tes}
    order = ['tss', 'ic', 'tes']

    c_dict, order = rm_color_cats(c_dict, order, cats)
    return c_dict, order

def get_end_colors():
    c_dict, order = get_sector_colors(['tes', 'tss'])
    return c_dict, order

def get_edge_colors():
    """
    Get colors for introns and exons
    """
    c_dict = {'intron': '#CC79A7', 'exon': '#009E73'}
    return c_dict

def get_biosample_colors():
    """
    Get colors for each biosample
    """
    d = os.path.dirname(__file__)
    fname = '{}/../refs/biosample_colors.tsv'.format(d)
    df = pd.read_csv(fname, sep='\t')

    c_dict = {}
    for ind, entry in df.iterrows():
        c_dict[entry.biosample] = entry.color
    order = df.biosample.tolist()

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


def get_ad_colors():
    c_dict = {'healthy': '#bb8f8f',
              'AD': '#b5bd61',
              np.nan: '#000000'}
    order = ('healthy', 'AD', np.nan)
    return c_dict, order

def plot_det_gene_len(df, opref='figures/'):
    pass

def plot_cell_line_tissue_det_venn(df, how='gene',
                                   nov='Known',
                                   opref='figures/'):
    """
    Plot a venn diagram showing how many genes / transcripts
    are detected between the cell line and tissue datasets

    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
            or filtered (for gene)
        how (str): Either 'gene' or 'iso'
        nov (str): Novelty, either 'Known', 'NIC', or 'NNC'
        opref (str): Where to save output figures
    """
    sns.set_context('paper', font_scale=1.8)

    if how == 'gene':
        df = df.loc[df.gene_novelty == 'Known']
        dataset_cols = get_dataset_cols()
        df = df[['annot_gene_id']+dataset_cols]
        df = df.groupby('annot_gene_id').sum().reset_index()

    tissues = get_sample_datasets('tissue')
    cell_lines = get_sample_datasets('cell_line')

    df['tissue'] = df[tissues].sum(1).astype(bool)
    df['cell_line'] = df[cell_lines].sum(1).astype(bool)

    if how == 'iso':
        df = df[['transcript_novelty', 'tissue', 'cell_line']]
        df = df.loc[df.transcript_novelty == nov]
    else:
        print(df.head())
        df = df[['annot_gene_id', 'tissue', 'cell_line']]

    known_out_green = '#90D6C3'
    known_int_green = '#009E73'
    nnc_out_gold = '#F5DFAE'
    nnc_int_gold = '#E69F00'
    nic_out_orange = '#DEA67A'
    nic_int_orange = '#D55E00'

    if nov == 'Known':
        out_color = known_out_green
        int_color = known_int_green
    elif nov == 'NNC':
        out_color = nnc_out_gold
        int_color = nnc_int_gold
    elif nov == 'NIC':
        out_color = nic_out_orange
        int_color = nic_int_orange

    df = df.groupby(['tissue', 'cell_line']).count().reset_index()
    if how == 'iso':
        df.rename({'transcript_novelty': 'counts'}, axis=1, inplace=True)
    elif how == 'gene':
        df.rename({'annot_gene_id': 'counts'}, axis=1, inplace=True)
    print(df)
    intersection = df.loc[(df.cell_line)&(df.tissue), 'counts'].values[0]
    tissue = df.loc[(~df.cell_line)&(df.tissue), 'counts'].values[0]
    cell_line = df.loc[(df.cell_line)&(~df.tissue), 'counts'].values[0]

    counts = [cell_line, tissue, intersection]
    log_counts = [np.log2(n) for n in counts]
    log_counts = tuple(counts)

    v = venn2(subsets=log_counts, set_labels=('',''))
    v.get_patch_by_id('10').set_color(out_color)
    v.get_patch_by_id('01').set_color(out_color)
    v.get_patch_by_id('11').set_color(int_color)
    v.get_patch_by_id('10').set_edgecolor(int_color)
    v.get_patch_by_id('01').set_edgecolor(int_color)
    v.get_patch_by_id('11').set_edgecolor(int_color)
    v.get_patch_by_id('10').set_linewidth(5)
    v.get_patch_by_id('01').set_linewidth(5)
    v.get_patch_by_id('11').set_linewidth(5)
    v.get_patch_by_id('10').set_alpha(1)
    v.get_patch_by_id('01').set_alpha(1)
    v.get_patch_by_id('11').set_alpha(1)
    v.get_label_by_id('10').set_text(counts[0])
    v.get_label_by_id('01').set_text(counts[1])
    v.get_label_by_id('11').set_text(counts[2])

    fname = '{}{}_{}_venn.png'.format(opref, how, nov)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_region_widths(regions,
                       kind='annot',
                       opref='figures/human'):
    """
    Parameters:
        regions (dict of pandas DataFrame): Output from
            get_ic_tss_tes
        kind (str): Choose 'annot', 'obs', or 'all'
        opref (str): Output file prefix
    """

    # plot histogram of tss / tes region sizes
    sns.set_context('paper', font_scale=1.8)
    for c, temp in regions.items():
        ax = sns.displot(data=temp, x='len', kind='hist', linewidth=0, binwidth=50)
        ax.set(xlabel='{}s'.format(c.upper()))
        fname = '{}_{}_{}_region_widths.png'.format(opref, kind, c)
        plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_genes_n_ic_ends(counts,
                         kind='annot',
                         opref='figures/human'):
    """
    Parameters:
        counts (pandas DataFrame): DF output from get_ic_tss_tes
        kind (str): Choose from 'annot', 'all', 'obs'
        opref (str): Where to save
    """
    # plot # tss / tes vs ic, color by # genes
    counts['tss_tes'] = counts.tss+counts.tes
    counts.head()

    # calculate how many genes have n tss+tes and n ics
    temp = counts[['tss_tes', 'intron_chain', 'gid']].groupby(['tss_tes', 'intron_chain']).count().reset_index()
    temp.rename({'gid': 'n_genes'}, inplace=True, axis=1)
    temp['log10_n_genes'] = np.log10(temp.n_genes)

    # plot the figure
    sns.set_context('paper', font_scale=1.6)
    plt.figure(figsize=(6,8))

    # ax = sns.scatterplot(data=temp, x='tss_tes', y='intron_chain', hue='n_genes', size='n_genes', palette='viridis')
    ax = sns.scatterplot(data=temp, x='tss_tes', y='intron_chain', hue='log10_n_genes', size='log10_n_genes', palette='viridis')

    norm = plt.Normalize(temp['log10_n_genes'].min(), temp['log10_n_genes'].max())
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])

    # Remove the legend and add a colorbar
    ax.get_legend().remove()
    cb = ax.figure.colorbar(sm)
    cb.set_label('log10(# genes)')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = '# TSSs + # TESs'
    ylabel = '# intron chains'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel, xscale='log', yscale='log')

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    fname = '{}_{}_n_genes_ic_tss_tes.png'.format(opref, kind)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_ic_upset(h5,
                   subset=None,
                   sources=None,
                   gids=None,
                   max_n_subsets=None,
                   opref='figures/cerberus',
                   **kwargs):
    df, _, _, _, _, _ = read_h5(h5, as_pyranges=False)
    # filter ends for gene subset
    if subset:
        df = filter_cerberus_genes(df, subset=subset)

    if gids:
        df = df.loc[df.gene_id.isin(gids)]

    # get melted version of regions
    ic_upset = upsetplot.from_memberships(df.source.str.split(','), data=df)

    # filter for given sources
    if sources:
        temp = ic_upset.copy(deep=True)
        all_sources = temp.index.names
        temp = temp.reset_index()
        temp = temp.loc[temp[sources].any(axis=1)]
        drop_sources = list(set(all_sources)-set(sources))
        temp.drop(drop_sources, axis=1, inplace=True)
        temp.set_index(sources, inplace=True)
        ic_upset = temp.copy(deep=True)

    # make the plot
    c_dict = get_edge_colors()
    mode = 'intron'
    c = c_dict[mode]
    fig = plt.figure(figsize=(11,6))
    sns.set_context('paper', font_scale=1.5)
    upsetplot.plot(ic_upset, subset_size='auto',
                    show_counts='%d', sort_by='cardinality',
                    facecolor=c, fig=fig, shading_color='white', element_size=None,
                    **kwargs)

    fname = '{}_{}_source_upset.png'.format(opref, mode)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

    return ic_upset

def plot_end_upset(h5, mode,
                   subset=None,
                   sources=None,
                   gids=None,
                   max_n_subsets=None,
                   opref='figures/cerberus',
                   **kwargs):
    _, tss, tes, _, _, _ = read_h5(h5, as_pyranges=False)
    if mode == 'tss':
        df = tss
    elif mode == 'tes':
        df = tes

    # filter ends for gene subset
    if subset:
        df = filter_cerberus_genes(df, subset=subset)

    if gids:
        df = df.loc[df.gene_id.isin(gids)]

    # get melted version of regions
    end_upset = upsetplot.from_memberships(df.source.str.split(','), data=df)

    # limit just to n subsets
#     if max_n_subsets:
#         content = from_contents(data)
#         uniques, counts = np.unique(content.index, return_counts=True)

#         sorted_uniques = [x for _, x in sorted(zip(counts, uniques), reverse=True)]

    # filter for given sources
    if sources:
        temp = end_upset.copy(deep=True)
        all_sources = temp.index.names
        temp = temp.reset_index()
        temp = temp.loc[temp[sources].any(axis=1)]
        drop_sources = list(set(all_sources)-set(sources))
        temp.drop(drop_sources, axis=1, inplace=True)
        temp.set_index(sources, inplace=True)
        end_upset = temp.copy(deep=True)


    # make the plot
    c_dict, _ = get_end_colors()
    c = c_dict[mode]
    fig = plt.figure(figsize=(11,6))
    sns.set_context('paper', font_scale=1.5)
    upsetplot.plot(end_upset, subset_size='auto',
                    show_counts='%d', sort_by='cardinality',
                    facecolor=c, fig=fig, shading_color='white', element_size=None,
                   **kwargs)

    fname = '{}_{}_source_upset.png'.format(opref, mode)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

    return end_upset

def plot_n_ic_tss_tes(counts,
                      label_genes=None,
                      kind='annot',
                      opref='figures/human'):
    """
    Parameters:
        counts (pandas DataFrame): DF output from get_ic_tss_tes
        label_genes (list of str): List of gene names
        kind (str): Choose from 'annot', 'all', 'obs'
        opref (str): Where to save thing
    """

    xs = ['tss', 'tss', 'tes']
    ys = ['intron_chain', 'tes', 'intron_chain']
    hues = ['tes', 'intron_chain', 'tss']

    for x, y, hue in zip(xs, ys, hues):

        # plot the figure
        sns.set_context('paper', font_scale=1.6)
        plt.figure(figsize=(6,8))

        ax = sns.scatterplot(data=counts, x=x, y=y, hue=hue, s=20, palette='viridis')

        norm = plt.Normalize(counts[hue].min(), counts[hue].max())
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
        sm.set_array([])

        # Remove the legend and add a colorbar
        ax.get_legend().remove()
        cb = ax.figure.colorbar(sm)
        if hue == 'tss' or hue == 'tes':
            cb.set_label('# {}s'.format(hue.upper()))
        else:
            cb.set_label('# {}s'.format(hue))

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        if x == 'tss' or x == 'tes':
            xlabel = '# {}s'.format(x.upper())
        else:
            xlabel = '# {}s'.format(x)
        if y == 'tss' or y == 'tes':
            ylabel = '# {}s'.format(y.upper())
        else:
            ylabel = '# {}s'.format(y)

        # annotate genes that are kinda interesting
        if label_genes:
            xlim = ax.get_xlim()[1]
            ylim = ax.get_ylim()[1]
            for g in label_genes:
                if g in counts.gname.tolist():
                    x_txt = counts.loc[counts.gname == g, x].values[0]+(1/80)*xlim
                    y_txt = counts.loc[counts.gname == g, y].values[0]-(1/80)*ylim
                    plt.annotate(g, (x_txt,y_txt), fontsize='small', fontstyle='italic')
        _ = ax.set(xlabel=xlabel, ylabel=ylabel)

        fname = '{}_{}_{}_{}_{}_scatter.png'.format(opref, x,y,hue, kind)
        plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_n_reps_per_biosamp(df,
                            sample='cell_line',
                            opref='figures/'):
    """
    Plot a bar plot showing the number of libraries
        that went into each sample

    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
        sample (str): Either "tissue", "cell_line"
        opref (str): Output prefix to save figure
    """

    dataset_cols = get_sample_datasets(sample)
    df = df[dataset_cols]
    df = df.transpose()
    df.reset_index(inplace=True)
    df.rename({'index': 'dataset'}, axis=1, inplace=True)
    df['celltype'] = df.dataset.str.rsplit('_', n=2, expand=True)[0]

    if sample == 'tissue':
        # add in the tissue metadata
        d = os.path.dirname(__file__)
        fname = '{}/../refs/tissue_metadata.csv'.format(d)
    #     fname = '../../refs/tissue_metadata.csv'.format(d)
        tissue = pd.read_csv(fname)
        df = df.merge(tissue[['biosample', 'tissue']],
                        how='left', left_on='celltype',
                        right_on='biosample')
        df.drop('celltype', axis=1, inplace=True)
        df.rename({'tissue': 'celltype'}, axis=1, inplace=True)
        print('Found {} distinct tissues'.format(len(df.celltype.unique())))
    elif sample == 'cell_line':
        print('Found {} distinct cell lines'.format(len(df.celltype.unique())))

    df = df[['dataset', 'celltype']]

    ax = sns.countplot(data=df, x='celltype')
    ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if sample == 'tissue':
        ylabel = '# libraries'
        xlabel = 'Tissue'
    elif sample == 'cell_line':
        ylabel = '# libraries'
        xlabel = 'Cell line'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    fname = '{}{}_libs_per_biosamp.png'.format(opref, sample)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_major_iso_pis(sg, groupby,
                       opref='figures/human'):
    """
    Plots a histogram of pi values for the 1st-4th highest-expressed isoform
    per gene per sample.

    Parameters:
        sg (swan_vis SwanGraph): SwanGraph with abundance information
        groupby (str): Column in sg.adata.obs on which to calculat pi
        opref (str): Output file prefix
    """
    obs_col = groupby
    df, counts = swan.calc_pi(sg.adata, sg.t_df, obs_col=obs_col)
    df = df.transpose()

    # merge with gene info
    df = df.merge(sg.t_df['gid'],
                  how='left', left_index=True, right_index=True)

    # get top isoform per gene
    df.reset_index(inplace=True)
    df = df.melt(id_vars=['tid', 'gid'])

    df.rename({'variable': obs_col,
               'value': 'pi'}, axis=1, inplace=True)
    df = df.sort_values(by='pi', ascending=False)

    # remove unexpressed transcripts
    df = df.loc[df.pi != 0]

    # add # isoforms / gene / dataset to use for filtering later
    temp = df[['gid', obs_col, 'tid']].groupby(['gid', obs_col]).count().reset_index()
    temp.rename({'tid': 'iso_counts'}, axis=1, inplace=True)
    df = df.merge(temp, how='left', on=['gid', obs_col])


    for i in range(1,5):
        # https://stackoverflow.com/questions/36310564/pandas-second-max-value-per-group-in-dataframe
        temp = df.groupby(['gid', obs_col]).head(i).groupby(['gid', obs_col]).tail(1).copy(deep=True)

        # make sure we're only looking at isoforms that have at least
        # i expressed isoforms
        temp = temp.loc[temp.iso_counts >= i]

        # plot distribution of max pi value
        # per gene, per tissue / age combo
        sns.set_context('paper', font_scale=2)
        ax = sns.displot(data=temp, x='pi', linewidth=0)

        if i == 1:
            xlabel = 'Highest pi per gene per {}'.format(obs_col)
        elif i == 2:
            xlabel = '2nd highest pi per gene per {}'.format(obs_col)
        elif i == 3:
            xlabel = '3rd highest pi per gene per {}'.format(obs_col)
        elif i > 3:
            xlabel = '{}th highest pi per gene per {}'.format(i, obs_col)
        ylabel = 'Number of isoforms'
        ax.set(ylabel=ylabel, xlabel=xlabel, xlim=(0,100))



def plot_ranked_exon_counts(sg,
                          df,
                          gene,
                          min_tpm=1,
                          gene_subset='polya',
                          sample='all',
                          groupby='library',
                          nov=['Known', 'NIC', 'NNC'],
                          opref='figures/human'):

    """
    Plot the ranked counts per exon for a given gene subset according
    to the input expression threshold and novelty categories

    Parameters:
        sg (swan_vis SwanGraph): SwanGraph with data from corresponding
    """

    # determine which isoforms are actually detected
    df = get_det_table(df,
                       how='iso',
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       sample=sample,
                       groupby=groupby,
                       nov=nov)
    tids = df.columns.tolist()

    # get isoforms from target gene
    df = sg.t_df.loc[sg.t_df.gname == gene]
    df = swan.pivot_path_list(df, 'path')
    df = df.merge(sg.edge_df, how='left', left_on='edge_id', right_index=True, suffixes=(None, '_dupe'))
    df = df.merge(sg.t_df[['tname']], how='left', left_index=True, right_index=True)
    df.drop('edge_id_dupe', axis=1, inplace=True)
    df = df[['edge_id', 'edge_type']]
    df.reset_index(inplace=True)
    tids = list(set(tids)&set(df.tid.tolist()))

    # limit to only detected isoforms and their edges
    df = df.loc[df.tid.isin(tids)]
    eids = df.edge_id.astype('int').tolist()
    df = sg.get_edge_abundance(kind='counts')
    df = df.loc[df.edge_id.isin(eids)]

    # sum up over datasets
    df.set_index('edge_type', inplace=True)
    df = df[sg.datasets]
    df['total_counts'] = df.sum(1)
    df.drop(sg.datasets, axis=1, inplace=True)

    # rank according to exp
    df = df.sort_values(by='total_counts', ascending=True)
    df.reset_index(inplace=True)

    sns.set_context('paper', font_scale=1.8)
    c_dict = get_edge_colors()

    for e in df.edge_type.unique():
        temp = df.loc[df.edge_type == e]
        temp['rank'] = [i for i in range(len(temp.index))]
        c = c_dict[e]

        ax = sns.catplot(data=temp, x='rank', y='total_counts', kind='bar',
                         linewidth=0, saturation=1, color=c)

        xlabel = 'Ranked $\it{}$ {}s'.format(gene, e)
        ylabel = 'Total counts'

        _ = ax.set(xlabel=xlabel, ylabel=ylabel, xticks=[])
        fname = '{}_{}_counts_per_{}.png'.format(opref, gene, e)
        plt.savefig(fname, dpi=300, bbox_inches='tight')

    return df

def plot_exon_hist(sg,
                   df,
                   gene,
                   min_tpm=1,
                   gene_subset='polya',
                   sample='all',
                   groupby='library',
                   nov=['Known', 'NIC', 'NNC'],
                   opref='figures/human'):
    """
    Plot a histogram of introns and exons and how many isoforms
    they're used in in a particular gene.

    Parameters:
        sg (swan_vis SwanGraph): SwanGraph of data
        df (pandas DataFrame): TALON abundance
    """
    # determine which isoforms are actually detected
    df = get_det_table(df,
                       how='iso',
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       sample=sample,
                       groupby=groupby,
                       nov=nov)
    tids = df.columns.tolist()

    # get isoforms from target gene
    df = sg.t_df.loc[sg.t_df.gname == gene]
    df = swan.pivot_path_list(df, 'path')
    df = df.merge(sg.edge_df, how='left', left_on='edge_id', right_index=True, suffixes=(None, '_dupe'))
    df = df.merge(sg.t_df[['tname']], how='left', left_index=True, right_index=True)
    df.drop('edge_id_dupe', axis=1, inplace=True)
    df = df[['edge_id', 'edge_type']]
    df.reset_index(inplace=True)
    tids = list(set(tids)&set(df.tid.tolist()))

    # limit to only detected isoforms
    df = df.loc[df.tid.isin(tids)]

    print('Found {} isoforms for {}'.format(len(tids), gene))

    c_dict = get_edge_colors()

    # groupby and count the number of isoforms that use each edge
    df = df.groupby(['edge_id', 'edge_type']).count().reset_index()
    df.rename({'tid':'n_isos'}, axis=1, inplace=True)

    sns.set_context('paper', font_scale=1.8)

    for e_type in df.edge_type.unique():
        temp = df.loc[df.edge_type == e_type]
        ax = sns.displot(data=temp, x='n_isos',
                         kind='hist',
                         color=c_dict[e_type], binwidth=5,
                         linewidth=0)

        xlabel = '# $\it{}$ isoforms'.format(gene)
        ylabel = '# {}s'.format(e_type)

        _ = ax.set(xlabel=xlabel, ylabel=ylabel)
        fname = '{}_{}_isos_per_{}.png'.format(opref, gene, e_type)
        plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_exons_per_iso(df,
                       nov=['Known'],
                       min_tpm=1,
                       gene_subset=None,
                       opref='figures/'):
    """
    Plots a boxplot showing # exons per isoform.

    Parameters:
        df (pandas DataFrame): TALON abundance file
        nov (list of str): Novelty categories to include
        min_tpm (float): Mininmum TPM to include isoform
        gene_subset (str): Choose from 'polya' or None
        opref (str): Output file prefix
    """

    t_df = df.copy(deep=True)

    # filter based on TPM, novelty, gene subset
    df, tids = get_tpm_table(df,
                             how='iso',
                             nov=nov,
                             min_tpm=min_tpm,
                             gene_subset=gene_subset)

    # add back in novelty and n exon info
    t_df = t_df[['annot_transcript_id', 'n_exons', 'transcript_novelty']]
    df = df.merge(t_df, how='left', left_index=True, right_on='annot_transcript_id')

    # plot the plot
    sns.set_context('paper', font_scale=1.6)
    plt.figure(figsize=(4,6))

    c_dict, order = get_talon_nov_colors(cats=nov)
    ax = sns.boxplot(data=df, x='transcript_novelty', y='n_exons',
                     order=order, palette=c_dict,
                     saturation=1)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = 'Transcript novelty'
    ylabel = '# exons'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    fname = '{}_exons_per_iso.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_gene_v_iso_sample_det(df,
                               sample='cell_line',
                               opref='figures/'):
    """
    Plot a hexbin density plot of the number of samples each gene and transcript is detected per librarygenes detected per library, by novelty type

    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
        sample (str): Either "tissue", "cell_line"
        opref (str): Output prefix to save figure
    """
    df = rm_sirv_ercc(df)
    dataset_cols = get_sample_datasets(sample)
    gridsize = int(len(dataset_cols)/2)
    print('Gridsize: {}'.format(gridsize))

    t_df = df[dataset_cols+['annot_transcript_id']].copy(deep=True)
    t_df.set_index('annot_transcript_id', inplace=True)
    t_df = t_df.astype(bool)
    t_df['n_samples_transcript'] = t_df.sum(1)
    t_df = t_df['n_samples_transcript'].to_frame()
    t_df.reset_index(inplace=True)

    g_df = df.loc[df.gene_novelty == 'Known'].copy(deep=True)
    g_df = g_df[dataset_cols+['annot_gene_id']]
    g_df = g_df.groupby('annot_gene_id').sum()
    g_df = g_df.astype(bool)
    g_df['n_samples_gene'] = g_df.sum(1)
    g_df = g_df['n_samples_gene'].to_frame()
    g_df.reset_index(inplace=True)

    cols = ['annot_gene_id', 'annot_transcript_id', 'transcript_novelty']
    t_df = t_df.merge(df[cols], how='left', on='annot_transcript_id')
    t_df = t_df.merge(g_df, how='left', on='annot_gene_id')

    # c_dict, order = get_talon_nov_colors(['Known', 'NIC', 'NNC'])
    c_dict, order = get_talon_nov_colors()
    sns.set_context('paper', font_scale=1.6)

    for nov in ['Known', 'NIC', 'NNC']:
        temp = t_df.loc[t_df.transcript_novelty == nov].copy(deep=True)
        ax = sns.jointplot(data=temp, x='n_samples_transcript',
                         y='n_samples_gene',
                         kind='hex',
                         color=c_dict[nov],
                         bins='log',
                         gridsize=gridsize)
        ax = ax.ax_joint

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        if sample == 'tissue':
            ylabel = '# tissue libraries gene was detected in'
            xlabel = '# tissue libraries transcript was detected in'
        elif sample == 'cell_line':
            ylabel = '# cell line libraries gene was detected in'
            xlabel = '# cell line libraries transcript was detected in'

        _ = ax.set(xlabel=xlabel, ylabel=ylabel)

        fname = '{}{}_gene_v_{}_iso_n_samp_det.png'.format(opref, sample, nov)
        plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_n_libs_v_avg_isos(df,
                           color='blue',
                           min_tpm=1,
                           gene_subset='polya',
                           sample='all',
                           nov=['Known', 'NIC', 'NNC'],
                           opref='figures/'):
    """
    Plot a scatterplot with a regression line for the average
        number of isos / gene vs. # libraries / sample

    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
        filt_df (pandas DataFrame): TALON abundance, filtered
        color (str): Color to plot plot in
        min_tpm (float): Minimum TPM to call a gene / iso as detected
        gene_subset (str): Subset of genes to use, 'polya' or None
        sample (str): Either "tissue", "cell_line", or None
        groupby (str): Either "sample", or "library",
            used to groupby datasets displayed
        nov (str): Novelty category of
            isoforms to consider
        opref (str): Output prefix to save figure
    """

    # get number of libraries per sample
    n_libs = get_n_libs_per_sample()

    # get avg isos
    df = get_isos_per_gene(df,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       sample=sample,
                       groupby='sample',
                       nov=nov)
    df = df.mean().to_frame().rename({0: 'avg_isos'}, axis=1)
    df.reset_index(inplace=True)
    df.rename({'index': 'biosample'}, axis=1, inplace=True)

    # merge with number of libraries
    df = df.merge(n_libs, how='left', on='biosample')

    sns.set_context('paper', font_scale=1.6)
    ax = sns.scatterplot(data=df, x='n_libraries', y='avg_isos', color='b')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = 'Number of libraries / sample'
    ylabel = 'Average # isoforms / gene'

    # get coeffs of linear fit
    c1 = 'n_libraries'
    c2 = 'avg_isos'
    slope, intercept, r_value, p_value, std_err = stats.linregress(df[c1],df[c2])
    lines = mpl.lines.Line2D([0], [0])
    label = 'm={0:.1f}'.format(slope)

    r2, pval_spear = stats.spearmanr(df[c1],df[c2])

    print('Slope of correlation: {}'.format(slope))
    print('R of correlation: {}'.format(r_value))
    print('R2 of correlation: {}'.format(r2))

    sns.regplot(data=df, x=c1, y=c2,
                scatter=False, ax=ax,
                color='b')
    sns.regplot(data=df, x=c1, y=c2,
        scatter=False, ax=ax, ci=0, color='b',
        line_kws={'linestyle':'-',
                  'label':"m={0:.1f}".format(slope)})

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    fname = '{}_libs_v_avg_isos.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_n_isos_gene_sample_dist():
    """
    Plot the distribution of # isoforms detected / gene / sample
    """

def plot_n_reads_v_avg_isos(df,
                            filt_df,
                            color='blue',
                            min_tpm=1,
                            gene_subset='polya',
                            sample='all',
                            groupby='sample',
                            nov=['Known', 'NIC', 'NNC'],
                            opref='figures/'):
    """
    Plot a scatterplot with a regression line for the average
        number of isos / gene vs. # reads / sample

    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
        filt_df (pandas DataFrame): TALON abundance, filtered
        color (str): Color to plot plot in
        min_tpm (float): Minimum TPM to call a gene / iso as detected
        gene_subset (str): Subset of genes to use, 'polya' or None
        sample (str): Either "tissue", "cell_line", or None
        groupby (str): Either "sample", or "library",
            used to groupby datasets displayed
        nov (str): Novelty category of
            isoforms to consider
        opref (str): Output prefix to save figure
    """

    # get number of reads from unfiltered data
    reads = get_reads_per_sample(df, groupby='sample')

    # get avg isos
    df = get_isos_per_gene(filt_df,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       sample=sample,
                       groupby=groupby,
                       nov=nov)
    df = df.mean().to_frame().rename({0: 'avg_isos'}, axis=1)
    df.reset_index(inplace=True)
    df.rename({'index': 'biosample'}, axis=1, inplace=True)

    # merge with read depth
    df = df.merge(reads, how='left', on='biosample')
    df.head()

    sns.set_context('paper', font_scale=1.6)
    ax = sns.scatterplot(data=df, x='n_reads', y='avg_isos', color='b')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = 'Reads per cell line / tissue'
    ylabel = 'Average # isoforms / gene'

    # get coeffs of linear fit
    c1 = 'n_reads'
    c2 = 'avg_isos'
    slope, intercept, r_value, p_value, std_err = stats.linregress(df[c1],df[c2])
    lines = mpl.lines.Line2D([0], [0])
    label = 'm={0:.1f}'.format(slope)

    r2, pval_spear = stats.spearmanr(df[c1],df[c2])

    print('Slope of correlation: {}'.format(slope))
    print('R of correlation: {}'.format(r_value))
    print('R2 of correlation: {}'.format(r2))

    sns.regplot(data=df, x=c1, y=c2,
                scatter=False, ax=ax,
                color='b')
    sns.regplot(data=df, x=c1, y=c2,
        scatter=False, ax=ax, ci=0, color='b',
        line_kws={'linestyle':'-',
                  'label':"m={0:.1f}".format(slope)})

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    fname = '{}_reads_v_avg_isos.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_gene_v_iso_det(df, filt_df,
                        sample='cell_line',
                        opref='figures/'):
    """
    Plot a scatterplot of the number of genes detected per library
        versus the number of transcripts detected, by novelty type

    Parameters:
        df (pandas DataFrame): TALON abundance, unfiltered
        filt_df (pandas DataFrame): TALON abundance, filtered
        sample (str): Either "tissue", "cell_line"
        opref (str): Output prefix to save figure
    """

    df = rm_sirv_ercc(df)
    filt_df = rm_sirv_ercc(filt_df)
    dataset_cols = get_sample_datasets(sample)

    # only known genes
    gene_df = df.loc[df.gene_novelty == 'Known'].copy(deep=True)
    gene_df = gene_df[dataset_cols+['annot_gene_id']]
    gene_df = gene_df.groupby('annot_gene_id').sum()

    gene_df = gene_df.astype(bool)
    ind_cols = gene_df.columns
    gene_df = gene_df.transpose()
    gene_df['n_genes'] = gene_df.sum(1)
    gene_df = gene_df['n_genes'].to_frame()

    # known transcripts
    t_df = filt_df.loc[filt_df.transcript_novelty == 'Known'].copy(deep=True)
    t_df = t_df[dataset_cols+['annot_transcript_id']]
    t_df = t_df.astype(bool)
    ind_cols = t_df.columns
    t_df = t_df.transpose()
    t_df['Known'] = t_df.sum(1)
    t_df = t_df['Known'].to_frame()

    gene_df = gene_df.merge(t_df, left_index=True, right_index=True)

    # nic transcripts
    t_df = filt_df.loc[filt_df.transcript_novelty == 'NIC'].copy(deep=True)
    t_df = t_df[dataset_cols+['annot_transcript_id']]
    t_df = t_df.astype(bool)
    ind_cols = t_df.columns
    t_df = t_df.transpose()
    t_df['NIC'] = t_df.sum(1)
    t_df = t_df['NIC'].to_frame()

    gene_df = gene_df.merge(t_df, left_index=True, right_index=True)

    # nnc transcripts
    t_df = filt_df.loc[filt_df.transcript_novelty == 'NNC'].copy(deep=True)
    t_df = t_df[dataset_cols+['annot_transcript_id']]
    t_df = t_df.astype(bool)
    ind_cols = t_df.columns
    t_df = t_df.transpose()
    t_df['NNC'] = t_df.sum(1)
    t_df = t_df['NNC'].to_frame()

    gene_df = gene_df.merge(t_df, left_index=True, right_index=True)
    gene_df.reset_index(inplace=True)

    df = gene_df.melt(id_vars=['index', 'n_genes'],
                  value_vars=['Known','NIC','NNC'])
    df.rename({'value': 'transcript_counts'}, axis=1, inplace=True)
    df.rename({'variable': 'novelty'}, axis=1, inplace=True)

    c_dict, order = get_talon_nov_colors(['Known', 'NIC', 'NNC'])
    sns.set_context('paper', font_scale=1.6)

    ax = sns.jointplot(data=df, x='transcript_counts', y='n_genes',
                     hue='novelty', palette=c_dict,
    #                  xlim=(0,xlim), ylim=(0,ylim),
                     joint_kws={'data':df, 's':40, 'alpha':1})
    ax = ax.ax_joint

    ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_legend().remove()

    if sample == 'tissue':
        ylabel = 'Known genes per tissue library'
        xlabel = 'Transcripts per tissue library'
    elif sample == 'cell_line':
        ylabel = 'Known genes per cell line library'
        xlabel = 'Transcripts per cell line library'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    fname = '{}{}_gene_v_iso_det.png'.format(opref, sample)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_avg_isos_per_gene(df,
                           min_tpm=1,
                           gene_subset='polya',
                           sample='all',
                           groupby='sample',
                           nov=['Known', 'NIC', 'NNC'],
                           opref='figures/'):
    """
    Plot the average number of isoforms per gene that are seen in each
    sample or library.

    Parameters:
        df (pandas DataFrame): TALON abundance
        min_tpm (float): Minimum TPM to call a gene / iso as detected
        gene_subset (str): Subset of genes to use, 'polya' or None
        sample (str): Either "tissue", "cell_line", or None
        groupby (str): Either "sample", or "library",
            used to groupby datasets displayed
        nov (str): Novelty category of
            isoforms to consider
        opref (str): Output prefix to save figure

    Returns:
        df (pandas DataFrame): DataFrame detailing how many samples
            each gene or isoform was seen in
    """

    sns.set_context('paper', font_scale=1.6)

    # get # isos / gene / sample or library
    df = get_isos_per_gene(df,
                           min_tpm=min_tpm,
                           gene_subset=gene_subset,
                           sample=sample,
                           groupby=groupby,
                           nov=nov)

    # fill 0s (which represent unexpresed genes) with NaNs before
    # calculating the avg. # isos / sample
    df = df.replace(0, np.nan)

    # calculate the average to order the barplots
    avgs = df.mean().to_frame()
    avgs = avgs.sort_values(by=[0], ascending=False)
    order = avgs.index.tolist()

    # melt the df to get an entry for each sample and gene
    df = df.melt()
    df.rename({'variable': 'tissue',
               'value': 'n_isos'}, axis=1, inplace=True)
    df = df.loc[~df.n_isos.isnull()]

    # actually make the plot
    ax = sns.barplot(data=df, x='tissue', y='n_isos', order=order)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = 'Sample'
    ylabel = '# isoforms / gene'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.tick_params(axis="x", rotation=90, labelsize=10)


    fname = '{}_isos_per_gene_per_{}.png'.format(opref, sample)
    plt.savefig(fname, dpi=300, bbox_inches='tight')


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

def plot_corr(df, sample='cell_line',
              how='gene', nov='Known',
              opref='figures/', cluster=False):

    corrs = compute_corr(df, how=how, sample=sample)
    if cluster == False:
        g = sns.clustermap(corrs, cmap='viridis',
                        xticklabels=True, yticklabels=True,
                        row_cluster=False, col_cluster=False)
    else:
        g = sns.clustermap(corrs, cmap='viridis',
                xticklabels=True, yticklabels=True)

    if how == 'iso':
        if nov == 'Known':
            title_txt = 'Known transcript'
        else:
            title_txt = '{} transcript'.format(nov)
    elif how == 'gene':
        title_txt = 'Known gene'
    title = '{} expression correlations'.format(title_txt)

    g.fig.suptitle(title)

    fname = '{}{}_{}_correlation.pdf'.format(opref, nov, how)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

    fname = '{}{}_{}_correlation.png'.format(opref, nov, how)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

    fname = '{}{}_{}_correlation.tsv'.format(opref, nov, how)
    corrs.to_csv(fname, sep='\t')


def plot_ranked_biosamp(df, sample='cell_line', how='iso', nov='known',
                    ylim=None, opref='figures/'):
    sns.set_context('paper', font_scale=2)

    max_df, min_df = compute_detection(df, sample, how, nov)

    c_dict, order = get_talon_nov_colors()
    color = c_dict[nov]
    ax = sns.lineplot(data=max_df, x='rank',
                      y='n_cumulative', linewidth=3, color=color)
    ax = sns.lineplot(data=min_df, x='rank',
                      y='n_cumulative', linewidth=3,
                      ax=ax, color=color)

#     if how == 'gene':
#         # n gencode genes
#         ax.plot([0, 23], [58780,58780], color='gray',
#                 linestyle='--', linewidth=3)
#         font = {'color': 'gray', 'size': 16}
#         plt.text(9.5, 61500, "~60k GENCODE genes", fontdict=font)
#     if how == 'iso' and nov == 'Known':
#         # n gencode genes
#         ax.plot([0, 23], [206761,206761], color='gray',
#                 linestyle='--', linewidth=3)
#         font = {'color': 'gray', 'size': 16}
#         plt.text(5.75, 213000, "~200k GENCODE transcripts", fontdict=font)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if how == 'iso':
        if nov == 'Known':
            ylabel_txt = 'known isoforms'
        else:
            ylabel_txt = '{} isoforms'.format(nov)
    elif how == 'gene':
        ylabel_txt = 'known genes'

    ylabel = 'Cumulative {} detected'.format(ylabel_txt)
    xlabel = 'Ranked {}s'.format(sample)
    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    if ylim:
        _ = ax.set(ylim=(-100, ylim))

    fname = '{}cumulative_{}_{}_per_{}.pdf'.format(opref, nov, how, sample)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    fname = '{}cumulative_{}_{}_per_{}.png'.format(opref, nov, how, sample)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

    return max_df, min_df

def plot_short_long_det(df, opref, \
                    xlim, ylim, how='gene'):

    sns.set_context('paper', font_scale=2)

    if how == 'gene':
        c1 = 'n_genes'
        c2 = 'ill_gene_count'
    elif how == 'read':
        c1 = 'n_counts'
        c2 = 'ill_umi_count'

#     ax = sns.jointplot(data=df, x=c1, y=c2,
#                      xlim=(0,xlim), ylim=(0,ylim),
#                      joint_kws={'data':df, 's':40, 'alpha':1})
    ax = sns.jointplot(data=df, x=c1, y=c2,
                     xlim=(0,xlim), ylim=(0,ylim),
                     joint_kws={'data':df, 's':40, 'alpha':1})
    ax = ax.ax_joint

#     # plot regression lines and equation of regression lines
#     # https://stackoverflow.com/questions/48145924/different-colors-for-points-and-line-in-seaborn-regplot/68135585#68135585
#     # https://stackoverflow.com/questions/45902739/seaborn-annotate-the-linear-regression-equation
#     # https://stackoverflow.com/questions/62705904/add-entry-to-matplotlib-legend-without-plotting-an-object
#     lines = []
#     labels = []
#     for s in df['sample'].unique().tolist():
#         temp = df.loc[df['sample'] == s]
#         color = c_dict[s]
#         line_color = adjust_lightness(color, 0.5)

#         # get coeffs of linear fit
#         slope, intercept, r_value, p_value, std_err = stats.linregress(temp[c1],temp[c2])
#         lines += [mpl.lines.Line2D([0], [0], color=line_color)]
#         labels += ['m={0:.1f}'.format(slope)]

#         print('Slope of {} correlation: {}'.format(s, slope))

#         sns.regplot(data=temp, x=c1, y=c2,
#                     scatter=False, ax=ax, color=color)
#         sns.regplot(data=temp, x=c1, y=c2,
#             scatter=False, ax=ax, color=color, ci=0,
#             line_kws={'color':line_color,
#                       'linestyle':'-',
#                       'label':"m={0:.1f}".format(slope)})

    ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_legend().remove()

    if how == 'gene':
        _ = ax.set(xlabel='# genes/cell in PacBio', ylabel='# genes/cell detected in Illumina')
        plt.savefig('{}_genes_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')
    elif how == 'read':
        _ = ax.set(xlabel='# reads/cell in PacBio', ylabel='# UMIs/cell in Illumina')
    plt.savefig('{}_reads_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')


def plot_reads_per_bc(df, title, oprefix):
    """
    Parameters:
        df (pandas DataFrame): DataFrame of read_annot file
        title (str): Title of plot
        oprefix (str): Output file prefix
    """

    temp = get_reads_per_bc(df)

    sns.set_context('paper', font_scale=1.5)
    counts = temp['counts'].tolist()
    plt.plot(range(len(counts)),
            counts,
            color='lightgray',
            linewidth=2)
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_xlabel('Ranked cells by # reads (logscale)')
    ax.set_ylabel('# reads (logscale)')
    ax.set_title(title)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if title == 'Post-TALON':
        title = 'post_talon'

    plt.tight_layout()

    fname = '{}_{}_umis_v_barcodes.png'.format(oprefix, title)
    plt.savefig(fname)


def plot_gene_tpm_v_n_isos(df, filt_df,
                           min_tpm=1,
                           groupby='sample',
                           gene_subset='polya',
                           nov=['Known', 'NIC', 'NNC'],
                           opref='figures/human'):
    """
    Plot a scatterplot of gene tpm / library or sample vs.
    number of isoforms per gene / library or sample

    Parameters:
        df (pandas DataFrame): TALON abundance file
        filt_df (pandas DataFrame): filtered talon abundance file
    """

    tpm_df, _ = get_tpm_table(df,
                       how='gene',
                       min_tpm=min_tpm,
                       groupby=groupby,
                       gene_subset=gene_subset)

    iso_df = get_isos_per_gene(filt_df,
                               min_tpm=min_tpm,
                               gene_subset=gene_subset,
                               groupby=groupby,
                               nov=nov)
    tpm_df = tpm_df.melt(ignore_index=False, value_name='tpm')
    iso_df = iso_df.melt(ignore_index=False, value_name='n_iso', var_name='biosample').fillna(0)

    tpm_df.index.name = 'annot_gene_id'
    tpm_df.reset_index(inplace=True)
    iso_df.reset_index(inplace=True)

    df = tpm_df.merge(iso_df, how='outer', on=['annot_gene_id', 'biosample'])

    sns.set_context('paper', font_scale=1.6)
    ax = sns.scatterplot(data=df, x='n_iso', y='tpm')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = '# isoforms / gene / sample'
    ylabel = 'TPM / gene / sample'

    # _ = ax.set(xlabel=xlabel, ylabel=ylabel, yscale='log')
    _ = ax.set(xlabel=xlabel, ylabel=ylabel)


    fname = '{}_isos_v_tpm.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

    # add gname
    gene_df, _, _ = get_gtf_info(how='gene')
    gene_df = gene_df[['gid', 'gname']]
    df = df.merge(gene_df, how='left', left_on='annot_gene_id', right_on='gid')

    return df

def plot_det_vs_gencode_isos(df,
                         min_tpm=1,
                         gene_subset='polya',
                         nov=['Known', 'NIC', 'NNC'],
                         label_genes=None,
                         opref='figures/',
                         ver='v29',
                         ylim=None,
                         xlim=None):
    """
    Plot a scatterplot of the the total number of isoforms detected
    vs the number of isoforms in gencode
    """

    # detected isoforms
    det_df = get_isos_per_gene(df,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       groupby='all',
                       nov=nov)
    det_df.rename({'all': 'n_isos_det'}, axis=1, inplace=True)
    det_df.reset_index(inplace=True)
    det_df['gid'] = cerberus.get_stable_gid(det_df, 'annot_gene_id')

    # annotated isoforms
    gc_df = get_n_gencode_isos(subset='polya', ver=ver)
    gc_df.drop('gid', axis=1, inplace=True)
    gc_df.rename({'gid_stable': 'gid'}, axis=1, inplace=True)
    gc_df = gc_df[['gid', 'n_isos_gencode']]
    # pdb.set_trace()


    df = det_df.merge(gc_df, how='left', on='gid')

    # add gene name
    gene_df, _, _ = get_gtf_info(how='gene', subset='polya', ver=ver, add_stable_gid=True)
    gene_df.drop('gid', axis=1, inplace=True)
    gene_df.rename({'gid_stable': 'gid'}, axis=1, inplace=True)
    gene_df = gene_df[['gid', 'gname']]
    df = df.merge(gene_df, how='left', on='gid')

    # add a pseudocount of 1 to each metric
    df.n_isos_det = df.n_isos_det+1
    df.n_isos_gencode = df.n_isos_gencode+1

    # plot the figure
    sns.set_context('paper', font_scale=1.6)
    plt.figure(figsize=(6,6))
    ax = sns.scatterplot(data=df, x='n_isos_det', y='n_isos_gencode')

    fig = plt.gcf()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = 'Total # isoforms / gene'
    ylabel = '# isoforms / gene in GENCODE'
    _ = ax.set(xlabel=xlabel, ylabel=ylabel, xscale='log', yscale='log')

    fig = plt.gcf()

    # set x and y lims if provided
    if xlim:
        xlim = (0, xlim)
        ax.set(xlim=xlim)
    if ylim:
        ylim = (0, ylim)
        ax.set(ylim=ylim)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    print(lims)

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

#     fig = plt.gcf()
#     print('3')
#     print(fig.get_size_inches())

    # annotate genes that are kinda interesting
    if label_genes:
        xlim = ax.get_xlim()[1]
        ylim = ax.get_ylim()[1]
        for g in label_genes:
            if g in df.gname.tolist():
                # x = df.loc[df.gname == g, 'n_isos_det'].values[0]+math.log10((2/75)*xlim)
                # y = df.loc[df.gname == g, 'n_isos_gencode'].values[0]-math.log10((2/75)*ylim)
                x = df.loc[df.gname == g, 'n_isos_det'].values[0]
                y = df.loc[df.gname == g, 'n_isos_gencode'].values[0]
                if x > 0.2 and y > 0.2:
                    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.annotate.html#matplotlib.axes.Axes.annotate
                    plt.annotate(g, (x,y), fontsize='small', fontstyle='italic', xytext=(4,-5), textcoords='offset pixels')


#     fig = plt.gcf()
#     print('4')
#     print(fig.get_size_inches())

    fname = '{}_total_v_gencode_isos_per_gene.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

    # remove pseudocount for analysis
    df['n_isos_det'] = df.n_isos_det-1
    df['n_isos_gencode'] = df.n_isos_gencode-1

    return df

def plot_max_vs_all_isos(df,
                         min_tpm=1,
                         gene_subset='polya',
                         groupby='sample',
                         nov=['Known', 'NIC', 'NNC'],
                         label_genes=None,
                         opref='figures/',
                         ylim=None,
                         xlim=None):
    """
    Plot a scatterplot of the maximum number of isoforms detected
    per sample or library vs. the total number of isoforms detected
    """
    df_copy = df.copy(deep=True)

    # get maximum number of detected isoforms per library or sample
    max_isos = get_isos_per_gene(df,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       groupby=groupby,
                       nov=nov)
    max_isos = max_isos.max(axis=1).to_frame()
    max_isos.rename({0: 'max_isos'}, axis=1, inplace=True)

    # get total number of detected isoforms overall
    total = get_isos_per_gene(df_copy,
                       min_tpm=min_tpm,
                       gene_subset=gene_subset,
                       groupby='all',
                       nov=nov)
    total.rename({'all': 'total_isos'}, axis=1, inplace=True)

    # merge
    df = max_isos.merge(total, left_index=True, right_index=True)

    # add gene name
    gene_df, _, _ = get_gtf_info(how='gene')
    df = df.merge(gene_df, how='left', left_index=True, right_on='gid')

    # plot the figure
    sns.set_context('paper', font_scale=1.6)
    ax = sns.scatterplot(data=df, x='total_isos', y='max_isos')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # set x and y lims if provided
    if xlim:
        xlim = (0, xlim)
        ax.set(xlim=xlim)
    if ylim:
        ylim = (0, ylim)
        ax.set(ylim=ylim)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # annotate genes that are kinda interesting
    if label_genes:
        xlim = ax.get_xlim()[1]
        ylim = ax.get_ylim()[1]
        for g in label_genes:
            if g in df.gname.tolist():
                x = df.loc[df.gname == g, 'total_isos'].values[0]+(2/75)*xlim
                y = df.loc[df.gname == g, 'max_isos'].values[0]-(1.5/75)*ylim
                plt.annotate(g, (x,y), fontsize='small', fontstyle='italic')

    xlabel = 'Total # isoforms / gene'
    ylabel = 'Max. # isoforms / gene in one sample'
    _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    fname = '{}_max_v_all_isos_per_gene.png'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

    return df

def plot_isos_per_gene_hist(df,
                            min_tpm=1,
                            gene_subset='polya',
                            sample='all',
                            groupby='sample',
                            nov=['Known'],
                            rm_1=False,
                            pseudocount=None,
                            opref='figures/'):
    """
    Plots dist. of # isos / gene across the different samples

    Parameters:
        df (pandas DataFrame): TALON abundance file
        gene_subset (str): Choose from None or 'polya' or 'tf'
        sample (str): Choose from 'all', 'tissue', or 'cell_line'
        min_tpm (float): Min. TPM val for at least one library
        groupby (str): Choose from 'sample' or 'library'
        nov (list of str): Novelty categories to consider
        rm_1 (bool): Whether or not to remove 1-count gene / sample
        pseudocount (int): Add a pseudocount
        opref (str): Output file prefix
    """

    df = get_isos_per_gene(df,
                           min_tpm=min_tpm,
                           gene_subset=gene_subset,
                           groupby=groupby,
                           sample=sample,
                           nov=nov)

    # get long form dataframe
    df = df.melt(ignore_index=False)

    # remove 0-count gene / sample combos (as these
    # are genes that were not detected)
    df = df.loc[df.value != 0]

    # remove 1-count gene / sample combos too
    if rm_1:
        df = df.loc[df.value != 1]

    # add pseudocount if wanted
    if pseudocount:
        df.value = df.value + pseudocount

    sns.set_context('paper', font_scale=2)

    ax = sns.displot(data=df, x='value', kind='hist', binwidth=1, linewidth=0)
    xlabel = '# isoforms / gene / sample'
    ylabel = 'Number of genes'

    _ = ax.set(xlabel=xlabel, ylabel=ylabel, yscale='log')
    for a in ax.axes.flat:
        a.yaxis.set_major_formatter(ScalarFormatter())

    plt.savefig('{}_hist_isos_per_gene_per_sample.png'.format(opref), \
                dpi=300, bbox_inches='tight')
    return df, ax

def plot_det_len_kde(df,
                     how='gene',
                     subset='polya',
                     min_tpm=1,
                     xlim=None,
                     split_biotypes=False,
                     ver='v29',
                     opref='figures/'):
    """
    Plots dist. of gene or transcript length based on whether
    gene or transcript was detected at the given TPM value.

    Parameters:
        df (pandas DataFrame): TALON abundance file
        how (str): Choose from 'gene' or 'iso'
        subset (str): Choose from None or 'polya'
        min_tpm (float): Min. TPM val for at least one library
        xlim (float): Maximum length to display
        split_biotypes (bool): Split detected and undetected
            transcripts by biotype
        opref (str): Output file prefix

    Returns:
        df (pandas DataFrame): DataFrame used to plot from
    """
    df, ids = get_tpm_table(df,
                   how=how,
                   min_tpm=min_tpm,
                   nov=['Known'],
                   gene_subset=subset)
    gene_df, _, _ = get_gtf_info(how=how, subset=subset, ver=ver)
    df.reset_index(inplace=True)

    if how == 'gene':
        col = 'annot_gene_id'
        ref_col = 'gid'
        x_col = 'length'
    elif how == 'iso':
        col = 'annot_transcript_id'
        ref_col = 'tid'
        x_col = 't_len'
    df = df[col].to_frame()

    df = df.merge(gene_df, how='outer',
                  left_on=col, right_on=ref_col)

    df['detected'] = True
    df.loc[df[col].isnull(), 'detected'] = False

    sns.set_context('paper', font_scale=2)

    ax = sns.displot(data=df, x=x_col, kind='kde',
                     linewidth=3, hue='detected', common_norm=True)

    if how == 'gene':
        xlabel = 'Gene length'
        ylabel = 'Detected GENCODE genes'
    elif how == 'iso':
        xlabel = 'Transcript length'
        ylabel = 'Detected GENCODE transcripts'

    if xlim:
        _ = ax.set(xlabel=xlabel, ylabel=ylabel, xlim=(0,xlim))
    else:
        _ = ax.set(xlabel=xlabel, ylabel=ylabel)

    plt.savefig('{}_{}_det_{}_tpm_length.pdf'.format(opref, how, min_tpm), \
                dpi=300, bbox_inches='tight')
    return df

def plot_read_len_kde(df, opref):
    sns.set_context('paper', font_scale=2)

#     ax = sns.displot(data=df, x='read_length', hue=hue,
#                  palette=c_dict, kind='kde', hue_order=order, linewidth=3,
#                  common_norm=common_norm)
    ax = sns.displot(data=df, x='read_length', kind='kde', linewidth=3)
#     ax.set(xlabel='Read length', ylabel='KDE of reads',
#           title='Length distribution of Reads', xlim=(0,7500),
#           xticks=[0, 2500, 5000, 7500])
    plt.savefig('{}_read_length_kde.pdf'.format(opref), dpi=300, bbox_inches='tight')

def plot_cluster_proportions(adata,
                             bar_key,
                             color_key,
                             opref='figures/'):

    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([color_key, bar_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index()
    props = props.pivot(columns=bar_key, index=color_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)

    fig, ax = plt.subplots(dpi=300)
    fig.patch.set_facecolor("white")
    fig.set_size_inches(15, 5)

    cmap = None
    cluster_palette = '{}_colors'.format(color_key)
    if cluster_palette in adata.uns.keys():
        cluster_palette = adata.uns[cluster_palette]
        cmap = sns.palettes.blend_palette(
            cluster_palette,
            n_colors=len(cluster_palette),
            as_cmap=True)

    props.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        legend=None,
        colormap=cmap
    )

    ax.legend(bbox_to_anchor=(1.1, 2), title=color_key)
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=90)
    ax.set_xlabel(props.index.name.capitalize())
    ax.set_ylabel("Proportion")
    ax.grid(False)

    fname = opref+'{}_by_{}_prop.png'.format(bar_key, color_key)
    plt.savefig(fname, bbox_inches='tight')

def plot_depth_by_tech(adata, how, opref,
                       hue=None, xlim=None, ylim=None):
    df = adata.obs
    sns.set_context('paper', font_scale=2)

    if how == 'gene':
        c1 = 'n_genes'
        c2 = 'sr_gene_count'
    elif how == 'read':
        c1 = 'n_counts'
        c2 = 'sr_umi_count'

    if xlim:
        xlim = (0, xlim)
    if ylim:
        ylim = (0, ylim)
    ax = sns.jointplot(data=df, x=c1, y=c2,
                       xlim=xlim, ylim=ylim,
                       joint_kws={'data':df, 's':40, 'alpha':1},
                       hue=hue)
    ax = ax.ax_joint

    # ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # ax.get_legend().remove()

    if how == 'gene':
        _ = ax.set(xlabel='# genes/cell in PacBio', ylabel='# genes/cell detected in Illumina')
        plt.savefig('{}_genes_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')
    elif how == 'read':
        _ = ax.set(xlabel='# reads/cell in PacBio', ylabel='# UMIs/cell in Illumina')
    plt.savefig('{}_reads_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')

def add_perc(ax, data, feature):
    total = data[feature].sum()
    ylim = ax.get_ylim()[1]
    n_cats = len(ax.patches)
    for p in ax.patches:
        percentage = '{:.1f}%'.format(100 * p.get_height()/total)
#         x = p.get_x() + p.get_width() / 2 - 0.45
        x = p.get_x() + p.get_width() / 2 - (0.065)*n_cats
        y = p.get_y() + p.get_height() + ylim*0.00625
        ax.annotate(percentage, (x, y), size = 12)

def plot_read_novelty(df, opref, c_dict, order,
                      ylim=None, title=None,
                      datasets='all'):
    sns.set_context("paper", font_scale=1.6)

    temp = df.copy(deep=True)

    # filter on datasets
    if datasets != 'all':
        temp = temp.loc[temp.dataset.isin(datasets)]

    # count number of reads per cat
    temp = temp[['transcript_novelty', 'read_name']].groupby('transcript_novelty').count()
    temp.reset_index(inplace=True)
    temp.rename({'read_name':'counts'}, axis=1, inplace=True)
    print(temp)

    # actual plotting
    g = sns.catplot(data=temp, x='transcript_novelty',
                y='counts', kind='bar',
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]
    g.set_ylabels('Reads')
    g.set_xlabels('Transcript novelty')

    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, temp, 'counts')

    if ylim:
        g.set(ylim=(0,ylim))

    # add title
    if not title:
        g.fig.suptitle('Reads per novelty category')
    else:
        g.fig.suptitle('{} reads per novelty category'.format(title))

    # save figure
    fname = '{}_read_novelty'.format(opref)
    g.savefig(fname+'.pdf', dpi=300)


def plot_supported_feats(filt_ab,
                         h5,
                         feat,
                         obs_source,
                         ref_sources,
                         support_sources,
                         opref,
                         **kwargs):
    """
    Plot a bar plot showing which observed features are supported
    from external data sources in the cerberus annotation

    Parameters:
        filt_ab (str): Path fo filtered abundance file
        h5 (str): Path to cerberus annotation h5 object
        feat (str): {'tss', 'tes', 'ic'}
        obs_source (str): Source in cerberus annotation to
            consider the "observed" data
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
    if feat == 'ic':
        c_key = 'splicing'
    else:
        c_key = feat
    temp_c_dict, order = get_sector_colors()
    c = temp_c_dict[c_key]
    if support_sources:
        order = ['Known', 'Supported', 'Novel']
    else:
        order = ['Known', 'Novel']
    c_dict, order = get_shade_colors(c, order)

    # plotting
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    plt.figure(figsize=(3,4))

    ax = sns.barplot(data=df, y='counts', x='support',
                     palette=c_dict, order=order,
                     saturation=1)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    xlabel = ''
    ylabel = '# observed {}s'.format(feat.upper())

    _ = ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.tick_params(axis="x", rotation=45)

    labels = [item.get_text() for item in ax.get_xticklabels()]
    # labels[0] = 'Known'
    # labels[1] = 'Novel'
    ax.set_xticklabels(labels)


    fname = '{}_{}_support.png'.format(opref, feat)
    plt.savefig(fname, dpi=500, bbox_inches='tight')

    fname = '{}_{}_support.pdf'.format(opref, feat)
    plt.savefig(fname, dpi=500, bbox_inches='tight')


def plot_transcript_novelty(df, oprefix,
                            ylim=None,
                            title=None,
                            whitelist=None,
                            sample=None,
                            novs=None,
                            save_type='pdf'):
    """
    Plot number of transcripts per novelty category.

    Parameters:
        df (pandas DataFrame): TALON read annot file or
        oprefix (str): Place to save
        ylim (int): y limit of resultant plot
        title (str): Title of resultant plot
        whitelist (list of str): List of transcript IDs to retain
        sample (str): Choose from 'cell_line' or 'tissue'
        save_type (str): Choose from 'pdf' or 'png'
    """
    sns.set_context('paper', font_scale=1.6)

    temp = df.copy(deep=True)

    c_dict, order = get_talon_nov_colors(cats=novs)

    # remove transcripts that are not on whitelist
    if whitelist:
        beep = temp.loc[temp.transcript_ID.isin(whitelist)].index.tolist()
        if len(beep) == 0:
            beep = temp.loc[temp.annot_transcript_id.isin(whitelist)].index.tolist()
        temp = temp.loc[beep]

    # filter on datasets
    if sample:
        datasets = get_sample_datasets(sample)
    else:
        datasets = get_dataset_cols()
    cols = ['transcript_ID', 'transcript_novelty']
    temp = temp[cols+datasets]

    temp['total_counts'] = temp[datasets].sum(1)
    temp = temp.loc[temp.total_counts > 0]

    # count number of isoforms per cat
#     temp = temp[['transcript_ID', 'transcript_novelty', 'read_name']].groupby(['transcript_ID', 'transcript_novelty']).count()
#     temp.reset_index(inplace=True)
#     temp.drop('read_name', axis=1, inplace=True)
    temp = temp[['transcript_ID', 'transcript_novelty']]
    temp = temp.groupby('transcript_novelty').count()
    temp.reset_index(inplace=True)
    temp.rename({'transcript_ID': 'counts'}, axis=1, inplace=True)
    print(temp)
    novs = ['NIC', 'Known', 'ISM_rescue', 'NNC']
    complete = temp.loc[temp.transcript_novelty.isin(novs), 'counts'].sum(axis=0)
    print('Number of complete isoforms: {}'.format(complete))


    # actual plotting
    sns.set_context('paper', font_scale=1.8)
    plt.figure(figsize=(4,6))
    g = sns.catplot(data=temp, x='transcript_novelty',
                y='counts', kind='bar',
                saturation=1,
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]
    g.set_ylabels('Isoforms')
    g.set_xlabels('Transcript novelty')

    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, temp, 'counts')

    if ylim:
        g.set(ylim=(0,ylim))

#     # add title
#     if not title:
#         g.fig.suptitle('Transcript models per novelty category')
#     else:
#         g.fig.suptitle('{} transcript models per novelty category'.format(title))

    # save figure
    fname = '{}{}_isoform_novelty'.format(oprefix,sample)
    if save_type == 'png':
        g.savefig(fname+'.png', dpi=300, bbox_inches='tight')
    elif save_type == 'pdf':
        g.savefig(fname+'.pdf', dpi=300, bbox_inches='tight')

    plt.show()
    plt.clf()

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





def plot_transcript_novelty_per(df,
                                gene='ELN',
                                min_tpm=1,
                                gene_subset='polya',
                                groupby='sample',
                                nov=['Known', 'NIC', 'NNC'],
                                opref='figures/'):

    """
    Plot the # of detected isoforms / novelty category for a particular
    gene across samples or libraries
    """

    # get metadata about the transcripts that we need
    t_df = df.copy(deep=True)
    t_df = t_df[['annot_transcript_id', 'annot_gene_name',
                 'annot_gene_id', 'transcript_novelty']]

    df = get_det_table(df,
              how='iso',
              min_tpm=min_tpm,
              gene_subset=gene_subset,
              groupby=groupby,
              nov=nov)
    df = df.transpose()

    # isolate isoforms that belong to gene of interest
    df = df.merge(t_df, how='left', left_index=True, right_on='annot_transcript_id')
    df = df.loc[df.annot_gene_name == gene]

    # calc # isos detected per nov
    df.drop(['annot_gene_id', 'annot_gene_name'], axis=1, inplace=True)
    df = df.melt(id_vars=['annot_transcript_id', 'transcript_novelty'])
    df.rename({'variable': groupby,
               'value': 'detected'}, axis=1, inplace=True)
    df = df.groupby([groupby, 'transcript_novelty', 'detected']).count().reset_index()
    df.rename({'annot_transcript_id':'counts'}, axis=1, inplace=True)
    df = df.loc[df.detected == True]

    c_dict, order = get_talon_nov_colors(nov)

    # actual plotting
    sns.set_context('paper', font_scale=1.6)
    plt.figure(figsize=(4,6))
    g = sns.catplot(data=df,
                x=groupby,
                hue='transcript_novelty',
                y='counts', kind='bar',
                palette=c_dict, hue_order=order,
                saturation=1)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]

    g.set_ylabels('$\it{}$ isoforms'.format(gene))
    g.set_xlabels('Transcript novelty')

    fname = '{}_{}_novelty_per_{}.png'.format(opref, gene, groupby)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_transcript_novelty_per_1(df,
                                  gene='ELN',
                                  dataset='h9_chondro',
                                  min_tpm=1,
                                  gene_subset='polya',
                                  groupby='sample',
                                  nov=['Known', 'NIC', 'NNC'],
                                  opref='figures/'):

    # get metadata about the transcripts that we need
    t_df = df.copy(deep=True)
    t_df = t_df[['annot_transcript_id', 'annot_gene_name',
                 'annot_gene_id', 'transcript_novelty']]

    df = get_det_table(df,
              how='iso',
              min_tpm=min_tpm,
              gene_subset=gene_subset,
              groupby=groupby,
              nov=nov)

    # limit to only those detected in sample / library of interest
    df = df.loc[dataset].to_frame()
    df = df.loc[df[dataset] == True]

    # isolate isoforms that belong to gene of interest
    df = df.merge(t_df, how='left', left_index=True, right_on='annot_transcript_id')
    df = df.loc[df.annot_gene_name == gene]

    df = df[['annot_transcript_id', 'transcript_novelty']].groupby('transcript_novelty').count().reset_index()
    df.rename({'annot_transcript_id': 'counts'}, axis=1, inplace=True)
    c_dict, order = get_talon_nov_colors(nov)

    # actual plotting
    sns.set_context('paper', font_scale=1.8)
    plt.figure(figsize=(4,6))
    g = sns.catplot(data=df, x='transcript_novelty',
                y='counts', kind='bar',
                saturation=1,
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]
    g.set_ylabels('$\it{}$ isoforms in {}'.format(gene, dataset))
    g.set_xlabels('Transcript novelty')

    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, df, 'counts')

    fname = '{}_{}_{}_novelty.png'.format(opref, gene, dataset)
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    
def plot_brain_tissue_cell_line_umap(swan_file,
                                     filt_ab,
                                     min_tpm, 
                                     gene_subset,
                                     ofile):
    df = pd.read_csv(filt_ab, sep='\t')
    df, tids = get_tpm_table(df,
                   how='iso',
                   min_tpm=min_tpm,
                   gene_subset=gene_subset)

    sg = swan.read(swan_file)
    sg.gene_adata.var['gid_stable'] = cerberus.get_stable_gid(sg.gene_adata.var, col='gid')
    sg.gene_adata.var.reset_index(drop=True)
    sg.gene_adata.var.set_index('gid_stable', inplace=True)

    # transcript adata
    adata = sg.adata[:, tids]
    len(adata.var.index)

    # scanpy processing

    # normalize data matrix to 10,000 count
    sc.pp.normalize_total(adata, target_sum=1e4)

    # log 
    sc.pp.log1p(adata)
    adata.raw = adata

    # find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes = 70000, flavor = 'seurat',
        min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.scale(adata, max_value=10)

    # sc.tl.pca(adata, use_highly_variable=False)
    sc.tl.pca(adata, use_highly_variable=True)

    sc.pp.neighbors(adata, metric='cosine')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # add the brain_tissue_cell_line classification
    obs_col = 'brain_tissue_cell_line'
    adata.obs[obs_col] = adata.obs.classification.astype('str')
    adata.obs.loc[adata.obs['sample']=='brain', obs_col] = 'brain'

    # add colors
    # map values in order specific to
    cmap, order = get_tissue_cell_line_colors()
    cmap2, _ = get_biosample_colors()
    brain_color = cmap2['brain'] 
    order += ['brain']
    print(order)
    cmap['brain'] = brain_color
    print(cmap)
    adata.obs[obs_col] = adata.obs[obs_col].astype('category')
    obs_order = list(adata.obs_names)
    sample_order = adata.obs[obs_col].cat.categories.tolist()
    print(sample_order)
    print(adata.obs[obs_col].unique())
    sample_colors = [cmap[s] for s in sample_order]
    adata.uns['{}_colors'.format(obs_col)] = sample_colors

    # also store rgb values in dict for use with gen_report
    for key, item in cmap.items():
        item = item[1:]
        r,g,b = tuple(int(item[i:i+2], 16) for i in (0, 2, 4))
        cmap[key] = (r,g,b)
    adata.uns['{}_dict'.format(obs_col)] = cmap
    
    sc.set_figure_params(figsize=(3,3))

    sc.pl.umap(adata, color=obs_col, frameon=True, size=120, show=False)
    f = plt.gcf()
    f.savefig(ofile, dpi=500, bbox_inches='tight')

def plot_cell_line_tissue_read_len_v_ref(df, 
                                         gene_subset,
                                         ref_fname,
                                         xlim,
                                         ofile):
    
    cell_lines = get_sample_datasets('cell_line')
    tissues = get_sample_datasets('tissue')
    df['source'] = False
    df.loc[df.dataset.isin(cell_lines), 'source'] = 'Reads from cell lines'
    df.loc[df.dataset.isin(tissues), 'source'] = 'Reads from tissues'
    
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

def plot_sector_gene_counts(counts):
    temp = pd.DataFrame()
    for source in counts.source.unique():
        df = assign_gisx_sector(counts)
        df = df.loc[df.source == source]
        df = df[['gid', 'source', 'sector']].groupby(['source', 'sector']).count().reset_index()
        df.rename({'gid': 'n_genes'}, axis=1, inplace=True)
        df['total_genes'] = df.n_genes.sum()
        temp = pd.concat([temp, df])
    temp['perc'] = (temp.n_genes/temp.total_genes)*100
    temp = temp.loc[temp.sector != 'simple']

    y = '% of total genes'
    temp.rename({'perc': y}, axis=1, inplace=True)
    c_dict, order = get_sector_colors(['tss', 'splicing', 'tes', 'mixed', 'simple'])
    # plot both together
    sns.set_context('paper', font_scale=1.8)
    ax = sns.catplot(data=temp, x='source',
                y=y, hue='sector', kind='bar',
                palette=c_dict, saturation=1,
                hue_order=order)

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

    return temp


def plot_sankey(df,
                source,
                sink,
                counts,
                color,
                title):
    """
    Plot sankey diagram.

    Parameters:
        df (pandas DataFrame): DF w/ source, sink, and counts columns
        source (str): Column name for source from df
        sink (str): Column name for sink from df
        color (str): {'sector', 'nov'}
        title (str): Title for plot
    """

    if color == 'sector':
        c_dict, order = get_sector_colors()
    elif color == 'nov':
        c_dict, order = get_ic_nov_colors()
    print(c_dict)

    # df[source] = pd.Categorical(df[source], order)
    # df[sink] = pd.Categorical(df[sink], order)
    # df.sort_values([source, sink], inplace=True)

    # order.reverse()
    order_2 = order+order

    order_l = [o.capitalize() for o in order_2]
    order_l = [o if o != 'Tss' else 'TSS' for o in order_l]
    order_l = [o if o != 'Tes' else 'TES' for o in order_l]

    source_map = dict([(sect, i) for i, sect in enumerate(order)])
    sink_map = dict([(sect, i+len(order)) for i, sect in enumerate(order)])
    df['source'] = df[source].map(source_map)
    df['sink'] = df[sink].map(sink_map)

    def nodify(order):

        y_values = [i for i, beep in enumerate(order)]
        y_values += y_values
        x_values = [0 for i in range(len(order))]
        x_values += [1 for i in range(len(order))]

        y_values = [y/max(y_values) for y in y_values]
        x_values = [x/max(x_values) for x in x_values]
        x_values = [x if x > 0 else 0.01 for x in x_values]
        y_values = [y if y > 0 else 0.01 for y in y_values]

        return x_values, y_values

    ghost_cookie = nodify(order)

    nodes = dict(
        label=order_l,
        color=[c_dict[n] for n in order_2],
        x=ghost_cookie[0],
        y=ghost_cookie[1])

    # add coords to node labels
    # nodes['label'] = ['{}: ({},{})'.format(l, x, y) for l,x,y in zip(nodes['label'], nodes['x'], nodes['y'])]

    print(nodes)

    links = dict(
        source=df.source.tolist(),
        target=df.sink.tolist(),
        value=df[counts].tolist(),
        color=[c_dict[n] for n in df[source].tolist()]) # color links by source

    # print(links)

    data = go.Sankey(node=nodes, link=links, arrangement='snap')
    fig = go.Figure(data)
    fig.update_layout(title_text=title,
                      font_family='Times New Roman')
    fig.update_traces(textfont_family='Arial',
                      textfont_size=20,
                      selector=dict(type='sankey'))

    fig.show()
    return fig

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


def plot_browser_isos(ca, sg, gene, obs_col, obs_condition, filt_ab, major_set,
                              
                               h=0.1, w=56, x=14, fig_w=14, ref_sources=['v29','v40'], species='human',
                               h_space=None,
                               add_tss=False, add_ccre=False, major=False, order='expression',
                     light_shade=None, dark_shade=None):
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
        # ccre = ccre.intersect(region, strandedness=None, how='containment')
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

    def get_isos(ca, filt_ab, gene, sample):
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
        # get tpm df
        tpm_df = swan.calc_tpm(sg.adata, obs_col=obs_col).sparse.to_dense()
        tpm_df = tpm_df.transpose()
        tpm_df = tpm_df.loc[tids, obs_condition].to_frame()
        return tpm_df
    
    if major:
        tids = get_major_isos(major_set, gene, obs_condition)
    else:
        tids = get_isos(ca, filt_ab, gene, obs_condition)
    tpm_df = get_tpm_df(sg, tids, obs_col, obs_condition)

    # colormap definition
    # if cmap == 'gray':
    if not light_shade:
        light_shade = get_sector_colors()[0]['mixed']
    if not dark_shade:
        dark_shade = get_sector_colors()[0]['simple']
    cmap = mpl.colors.LinearSegmentedColormap.from_list('', [light_shade, dark_shade])
    # else:
    #     pass
        
    
    # font sizes
    # small_text = 6/(6.11/16)
    small_text = 20.3
    big_text = 6.71/(6.11/16)
    # print('small text size: {}'.format(small_text))
    # print('big text size: {}'.format(big_text))
    
    
    # height spacing b/w models
    if not h_space:
        h_space = h*1.75
    # print('h_space : {}'.format(h_space))
    # print('h: {}'.format(h))
    
    # plotting settings
    fig_len = len(tids)
    fig_len += 1 # for scale
    if add_tss: 
        fig_len += 1
    if add_ccre:
        fig_len += 1
    fig_h = h_space*(fig_len-1)+h
    fig_h += 0.3 # vertical spacing adjustment
    # print('fig h: {}'.format(fig_h))
    # print('fig w: {}'.format(fig_w))
    plt.figure(1, figsize=(fig_w, fig_h), frameon=False)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    ax = plt.gca()
    
    # plotting order
    tpm_df = tpm_df.sort_values(by=obs_condition, ascending=False)
    tpm_df = tpm_df[[obs_condition]]
    
    if order == 'tss':
        tpm_df = add_feat(tpm_df, kind='tss', col='index')
        tpm_df.sort_values(by=['tss', obs_condition], ascending=[True, False], inplace=True)
        tpm_df.drop('tss', axis=1, inplace=True)

    # x coords for gencode name and cerberus name
    x_gc = -4
    x_c = 11
    
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
    
    # label known transcripts with a * 
    tpm_df['Known'] = ''
    tpm_df = tpm_df.copy(deep=True)
    tpm_df = tpm_df.merge(ca.ic[['Name', 'novelty']], how='left', left_on='ic_id', right_on='Name')
    tpm_df.rename({'novelty':'ic_novelty'}, axis=1, inplace=True)
    tpm_df.drop('Name', axis=1, inplace=True)
    tpm_df = tpm_df.merge(ca.tss[['Name', 'novelty']], how='left', left_on='tss_id', right_on='Name')
    tpm_df.rename({'novelty':'tss_novelty'}, axis=1, inplace=True)
    tpm_df.drop('Name', axis=1, inplace=True)
    tpm_df = tpm_df.merge(ca.tes[['Name', 'novelty']], how='left', left_on='tes_id', right_on='Name')
    tpm_df.rename({'novelty':'tes_novelty'}, axis=1, inplace=True)
    tpm_df.drop('Name', axis=1, inplace=True)
    # tpm_df.loc[(tpm_df.tss_novelty=='Known')&\
    #            (tpm_df.ic_novelty=='Known')&\
    #            (tpm_df.tes_novelty=='Known'), 'Known'] = '*'
    if species == 'human':
        refs = ['v40', 'v29']
    elif species == 'mouse':
        refs = ['vM21', 'vM25']
    known_tids = ca.t_map.loc[ca.t_map.source.isin(refs)].transcript_id.unique().tolist()
    tpm_df.loc[tpm_df.transcript_id.isin(known_tids), 'Known'] = '*'
    
    # triplets rather than entire transcript name
    tpm_df['triplet'] = tpm_df.transcript_id.str.split('[', n=1, expand=True)[1]
    tpm_df['triplet'] = tpm_df.triplet.str.split(']', n=1, expand=True)[0]
    tpm_df['triplet'] = '['+tpm_df.triplet+']'
    
    i = 0
    # print('h: {}'.format(h))
    for index, entry in tpm_df.iterrows():
        
        # y coords
        y = (len(tpm_df.index) - i)*(h_space)
        y_text = y+(h/2)
        
        # tid
        tid = entry['transcript_id']
        gname = entry['gene_name']
        tname = entry['transcript_name']
        ref_tname = entry['original_transcript_name']
        trip = entry['triplet'] 
        iso_trip = '{} {}'.format(gname, trip)
        known = entry['Known']
    
        # color by TPM
        if len(tpm_df.index) == 1:
            norm_val = entry[obs_condition]
        else:
            norm_val = (entry[obs_condition]-tpm_df[obs_condition].min())/(tpm_df[obs_condition].max()-tpm_df[obs_condition].min())
        color = cmap(norm_val)
        ax = sg.plot_browser(tid, y=y, x=x, h=h, w=w, color=color, ax=ax) 
        # print('transcipt #{}, {}, y = {}'.format(i, iso_trip, y))

       # isoform name / novelty
        ax.text(x_c, y_text, known,
                verticalalignment='center', 
                horizontalalignment='center',
                size=small_text)   
        ax.text(x_gc, y_text, iso_trip,
                verticalalignment='center',
                horizontalalignment='left', 
                size=small_text) 
        
        i += 1
        
    # label the different columns
    # y = (len(tpm_df.index)+1)*(h_space)
   # # ax.text(x_c, y+(h/2), 'Cerberus Name',
    ##     verticalalignment='center', 
    ##     horizontalalignment='center')    
    ## ax.text(x_gc,y+(h/2), 'GENCODE v40 Name',
    ##         verticalalignment='center',
          #  # horizontalalignment='center')
    
    # ax.text(x_c, y+(h/2), 'Known',
    #     verticalalignment='center', 
    #     horizontalalignment='center', 
    #     size=big_text)    
    # ax.text(x_gc,y+(h/2), 'Isoform triplet',
    #         verticalalignment='center',
    #         horizontalalignment='center',
    #         size=big_text)   
    
    
    i = len(tpm_df.index)
    if add_tss:
        y = (len(tpm_df.index) - i)*(h_space)
        ax = plot_tss(ca, sg, tpm_df, x, y, h, ax)
        i += 1
    
    if add_ccre:
        y = (len(tpm_df.index) - i)*(h_space)
        ax = plot_ccre(ca, sg, x, y, h, ax)
        i += 1
        
    # add scale
    y = (len(tpm_df.index) - i)*(h_space)
    # print('y scale: {}'.format(y))
    ax = sg.pg.plot_scale(x, y, h, 14, ax)
        
    # remove axes
    plt.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False) 
    plt.tight_layout()
    
    # set x / y lim
    y_max = (len(tpm_df.index) - 0)*(h_space)+h
    y_min = y
    plt.ylim((y_min,y_max))
    plt.xlim((-4, 72))
    # print('ylim: ({},{})'.format(y_min,y_max))
    
    return ax, tpm_df

def plot_human_sample_legend(swan_file, ofile):
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    c_dict, order = get_biosample_colors()
    
    df = swan.read(swan_file).adata.obs.copy(deep=True)
    df = df[['sample', 'sample_display', 'classification']]
    df = df.drop_duplicates()
    df = df.sort_values(by=['classification', 'sample_display'], ascending=True)
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

    c_dict, order = get_lr_bulk_sample_colors()
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

def plot_mouse_brain_tissue_cell_line_umap(swan_file,
                                 filt_ab,
                                 min_tpm, 
                                 gene_subset,
                                 ofile):


    df = pd.read_csv(filt_ab, sep='\t')
    df, tids = get_tpm_table(df,
                   how='iso',
                   min_tpm=min_tpm,
                   species='mouse',
                   gene_subset=gene_subset)

    sg = swan.read(swan_file)
    adata = sg.adata[:, tids]


    # normalize data matrix to 10,000 count
    sc.pp.normalize_total(adata, target_sum=1e4)

    # log 
    sc.pp.log1p(adata)
    adata.raw = adata

    # find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes = 70000, flavor = 'seurat',
        min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.scale(adata, max_value=10)

    # sc.tl.pca(adata, use_highly_variable=False)
    sc.tl.pca(adata, use_highly_variable=True)

    sc.pp.neighbors(adata, metric='cosine')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    # make the figure square
    sc.set_figure_params(figsize=(3,3))

    # by cell line, tissue, brain
    adata = sg.adata[:, tids]
    len(adata.var.index)

    obs_col = 'brain_tissue_cell_line'
    adata.obs[obs_col] = 'tissue'
    adata.obs.loc[adata.obs['sample'].isin(['c2c12_myoblast', 'c2c12_myotube', 'f1219']), obs_col] = 'cell_line'
    adata.obs.loc[adata.obs['sample'].str.contains('hippo'), obs_col] = 'brain'
    adata.obs.loc[adata.obs['sample'].str.contains('cortex'), obs_col] = 'brain'

    # add colors
    # map values in order specific to
    cmap, order = get_tissue_cell_line_colors()
    cmap2, _ = get_biosample_colors()
    brain_color = cmap2['brain'] 
    order += ['brain']
    print(order)
    cmap['brain'] = brain_color
    print(cmap)
    adata.obs[obs_col] = adata.obs[obs_col].astype('category')
    obs_order = list(adata.obs_names)
    sample_order = adata.obs[obs_col].cat.categories.tolist()
    print(sample_order)
    print(adata.obs[obs_col].unique())
    sample_colors = [cmap[s] for s in sample_order]
    adata.uns['{}_colors'.format(obs_col)] = sample_colors

    # also store rgb values in dict for use with gen_report
    for key, item in cmap.items():
        item = item[1:]
        r,g,b = tuple(int(item[i:i+2], 16) for i in (0, 2, 4))
        cmap[key] = (r,g,b)
    adata.uns['{}_dict'.format(obs_col)] = cmap

    # normalize data matrix to 10,000 count
    sc.pp.normalize_total(adata, target_sum=1e4)

    # log 
    sc.pp.log1p(adata)
    adata.raw = adata

    # find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes = 70000, flavor = 'seurat',
        min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.scale(adata, max_value=10)

    # sc.tl.pca(adata, use_highly_variable=False)
    sc.tl.pca(adata, use_highly_variable=True)

    sc.pp.neighbors(adata, metric='cosine')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    sc.pl.umap(adata, color=obs_col, frameon=True, size=120, show=False)


    f = plt.gcf()
    f.savefig(ofile, dpi=500, bbox_inches='tight')
    
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
    
def plot_novel_supported_triplet_feats(filt_ab,
                                       h5,
                                       gene_subset,
                                       min_tpm,
                                       ofile):
                                       
    feats = ['tss', 'ic', 'tes']
    ref_sources = [['v29', 'v40'],
                   ['v29', 'v40'],
                   ['v29', 'v40']]
    support_sources = [['encode_cage', 'fantom_cage', 'encode_rampage'], 
                       ['gtex'],
                       ['pas', 'polya_atlas']]

    plt.figure(figsize=(3,20))
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    fig, axes = plt.subplots(1, len(feats), figsize=(10,4))
    i = 0

    for feat, ref_source, support_source in zip(feats, ref_sources, support_sources):
        print(feat)
        ax = axes[i]
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
                                   ofile):

    df = pd.read_csv(filt_ab, sep='\t')
    df, tids = get_tpm_table(df,
                   how='iso',
                   min_tpm=min_tpm,
                   gene_subset=gene_subset,
                   ic_nov=ic_nov,
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
                                             ofile):

    df = pd.read_csv(filt_ab, sep='\t')
    df, tids = get_tpm_table(df,
                   how='iso',
                   min_tpm=min_tpm,
                   gene_subset=gene_subset,
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


    
def plot_browser_isos_2(h5, 
                        swan_file, 
                        filt_ab,
                        major_isos,
                        gene, 
                        obs_col, 
                        obs_condition, 
                        ofile,
                        h=0.2,
                        ref_sources=['v29', 'v40'],
                        **kwargs):

    ca = cerberus.read(h5)
    sg = swan.read(swan_file)
    
    ax, tpm_df = plot_browser_isos(ca,
                                   sg,
                                   gene,
                                   obs_col,
                                   obs_condition,
                                   filt_ab,
                                   major_isos,
                                   h=h,
                                   ref_sources=ref_sources,
                                   **kwargs)
    plt.savefig(ofile, dpi=500)

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
    
    # get predominant transcripts in this sample
    prin_isos = pd.read_csv(major_isos, sep='\t')
    prin_isos = prin_isos.loc[prin_isos[obs_col] == obs_condition]
    prin_isos = prin_isos.loc[prin_isos.pi_rank == 1]
    prin_isos['gid_stable'] = cerberus.get_stable_gid(prin_isos, 'gid')
    prin_isos.drop(['gid', 'gname'], axis=1, inplace=True)
    prin_isos.rename({'gid_stable':'gid'}, axis=1, inplace=True)
    
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
    
# def plot_density_simplices(h5,
#                            sources,
#                            titles,
#                            gene_subset,
#                            ver):
    
#     def make_sector_source_bar_plots(counts, fname):
#         counts[['source', 'gid']].groupby('source').count()    
#         temp = pd.DataFrame()
#         for source in counts.source.unique():
#             df = assign_gisx_sector(counts)
#             df = df.loc[df.source == source]
#             df = df[['gid', 'source', 'sector']].groupby(['source', 'sector']).count().reset_index()
#             df.rename({'gid': 'n_genes'}, axis=1, inplace=True)
#             df['total_genes'] = df.n_genes.sum()
#             temp = pd.concat([temp, df])
#         temp['perc'] = (temp.n_genes/temp.total_genes)*100

#         cat1 = 'sector'
#         cat2 = 'source'
#         cat2_order = ['v40', 'obs_det', 'obs_major']
#         cat1_order = ['tss', 'splicing', 'tes', 'mixed', 'simple']
#         parent_c_dict, parent_order = get_sector_colors()

#         sns.set_context('paper', font_scale=2)
#         plt.figure(figsize=(2,20))
#         print('figsize changed')
#         fig, axes = plt.subplots(1, len(cat1_order), figsize=(20,3))

#         # loop through first category
#         for i, c1 in enumerate(cat1_order):
#             beep = temp.loc[temp[cat1] == c1]

#             c = parent_c_dict[c1]
#             c_dict, order = get_shade_colors(c, cat2_order)

#             # plotting
#             mpl.rcParams['font.family'] = 'Arial'
#             mpl.rcParams['pdf.fonttype'] = 42
#             # plt.figure(figsize=(3,4))

#             ax = sns.barplot(data=beep, y='perc', x=cat2,
#                              palette=c_dict, order=order,
#                              saturation=1, ax=axes[i])
#             ax.spines['right'].set_visible(False)
#             ax.spines['top'].set_visible(False)

#             xlabel = ''
#             if c1 in ['tss', 'tes']:
#                 c1 = c1.upper()
#             ylabel = '% genes in {} sector'.format(c1)

#             _ = ax.set(xlabel=xlabel, ylabel=ylabel)
#             ax.tick_params(axis="x", rotation=45)    
#             ax.set_xticklabels(['v40', 'Obs.', 'Obs. major'])

#             def add_perc_2(ax):
#                 ylim = ax.get_ylim()[1]
#                 n_cats = len(ax.patches)
#                 for p in ax.patches:
#                     percentage = '{:.1f}%'.format(p.get_height())
#                     x = p.get_x() + p.get_width() / 2
#                     y = p.get_y() + p.get_height() + ylim*0.01
#                     ax.annotate(percentage, (x, y), size=16, horizontalalignment='center')

#             add_perc_2(ax)
#         plt.savefig(fname, dpi=500, layout='tight', bbox_inches="tight")
        

#     # add biotype to subset on 
#     ca = cerberus.read(h5)
#     gene_df, _, _ = get_gtf_info(how='gene', ver=ver, add_stable_gid=True)
#     gene_df = gene_df[['gid_stable', 'biotype_category']]
#     gene_df.rename({'gid_stable': 'gid'}, axis=1, inplace=True)
#     ca.triplets = ca.triplets.merge(gene_df, how='left', on='gid')
    
#     gs_label = gene_subset
#     if gene_subset == 'polya':
#         gs = get_polya_cats()
#     else: 
#         gs = [gene_subset]

#     plot_df = pd.DataFrame()
#     for source, title in zip(sources, titles):

#         # subset dict
#         subset = {'source': source, 'biotype_category': gs}
#         print(subset)

#         # if we're looking at gencode, only take detected genes
#         if source == 'v40':
#             df = pd.read_csv(ab, sep='\t')
#             df, inds = get_tpm_table(df,
#                                      how='gene',
#                                      gene_subset=gs_label,
#                                      min_tpm=min_tpm)
#             subset['gid'] = inds

#         fname = 'figures/simplex_{}_{}.pdf'.format(source, gs_label)
#         temp = ca.plot_simplex(top='splicing_ratio',
#             subset=subset,
#             density=True,
#             density_scale=100,
#             density_cmap='Purples',
#             log_density=True,
#             sectors=True,
#             scatter=False,
#             legend=False,
#             title=title,
#             size_scale=0.2,
#             density_cbar=False,
#             fname=fname)

#         # add sectors from this subset to plot_df
#         temp = assign_gisx_sector(temp)
#         plot_df = pd.concat([temp, plot_df])

#     # create the bar plots
#     fname = 'figures/{}_genes_per_sector.pdf'.format(gs_label)
#     make_sector_source_bar_plots(plot_df, fname)
    
def plot_obs_human_simplex_with_centroid(h5, gene, fig_dir, **kwargs):
    ca = cerberus.read(h5)
    
    # add centroid as coord
    ca = get_centroids(ca,
                       **kwargs)
    
    # human
    c_dict, order = get_biosample_colors()
    c_dict[np.nan] = 'k'
    mmap = {'v40': '*', 'v29': 'x', 'obs_det': '^', 'cerberus': '2', 'sample_det_centroid': 'x'} 
    fname='{}/simplex_{}_centroid_vs_sample.pdf'.format(fig_dir, gene.lower())
    df = ca.plot_simplex(top='splicing_ratio', 
                gene=gene,
                hue='sample',
                cmap=c_dict,
                size='gene_tpm', 
                log_size=True,
                sectors=True,
                marker_style='source',
                mmap=mmap,
                legend=False,
                jitter=True,
                subset={'source': ['v40', 'obs_det', 'sample_det', 'sample_det_centroid']},
                size_scale=0.2,
                fname=fname)

def plot_density_simplices(h5,
                           ab,
                           min_tpm,
                           sources,
                           titles,
                           gene_subset,
                           ver,
                           species='human'):
    
    def make_sector_source_bar_plots(counts, fname, species='human'):
        counts[['source', 'gid']].groupby('source').count()    
        temp = pd.DataFrame()
        for source in counts.source.unique():
            df = assign_gisx_sector(counts)
            df = df.loc[df.source == source]
            df = df[['gid', 'source', 'sector']].groupby(['source', 'sector']).count().reset_index()
            df.rename({'gid': 'n_genes'}, axis=1, inplace=True)
            df['total_genes'] = df.n_genes.sum()
            temp = pd.concat([temp, df])
        temp['perc'] = (temp.n_genes/temp.total_genes)*100

        cat1 = 'sector'
        cat2 = 'source'
        if species == 'human':
            cat2_order = ['v40', 'obs_det', 'obs_major']
        elif species == 'mouse':
            cat2_order = ['vM25', 'obs_det', 'obs_major']
        cat1_order = ['tss', 'splicing', 'tes', 'mixed', 'simple']
        parent_c_dict, parent_order = get_sector_colors()

        sns.set_context('paper', font_scale=2)
        plt.figure(figsize=(2,20))
        print('figsize changed')
        fig, axes = plt.subplots(1, len(cat1_order), figsize=(20,3))

        # loop through first category
        for i, c1 in enumerate(cat1_order):
            beep = temp.loc[temp[cat1] == c1]

            c = parent_c_dict[c1]
            c_dict, order = get_shade_colors(c, cat2_order)

            # plotting
            mpl.rcParams['font.family'] = 'Arial'
            mpl.rcParams['pdf.fonttype'] = 42
            # plt.figure(figsize=(3,4))

            ax = sns.barplot(data=beep, y='perc', x=cat2,
                             palette=c_dict, order=order,
                             saturation=1, ax=axes[i])
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)

            xlabel = ''
            if c1 in ['tss', 'tes']:
                c1 = c1.upper()
            ylabel = '% genes in {} sector'.format(c1)

            _ = ax.set(xlabel=xlabel, ylabel=ylabel)
            ax.tick_params(axis="x", rotation=45)
            if species == 'human':
                ax.set_xticklabels(['v40', 'Obs.', 'Obs. major'])
            elif species == 'mouse':
                ax.set_xticklabels(['vM25', 'Obs.', 'Obs. major'])

            def add_perc_2(ax):
                ylim = ax.get_ylim()[1]
                n_cats = len(ax.patches)
                for p in ax.patches:
                    percentage = '{:.1f}%'.format(p.get_height())
                    x = p.get_x() + p.get_width() / 2
                    y = p.get_y() + p.get_height() + ylim*0.01
                    ax.annotate(percentage, (x, y), size=16, horizontalalignment='center')

            add_perc_2(ax)
        plt.savefig(fname, dpi=500, layout='tight', bbox_inches="tight")
        

    # add biotype to subset on 
    ca = cerberus.read(h5)
    gene_df, _, _ = get_gtf_info(how='gene', ver=ver, add_stable_gid=True)
    gene_df = gene_df[['gid_stable', 'biotype_category']]
    gene_df.rename({'gid_stable': 'gid'}, axis=1, inplace=True)
    ca.triplets = ca.triplets.merge(gene_df, how='left', on='gid')
    
    gs_label = gene_subset
    if gene_subset == 'polya':
        gs = get_polya_cats()
    else: 
        gs = [gene_subset]

    plot_df = pd.DataFrame()
    for source, title in zip(sources, titles):

        # subset dict
        subset = {'source': source, 'biotype_category': gs}
        print(subset)

        # if we're looking at gencode, only take detected genes
        if source in ['v40', 'vM25']:
            df = pd.read_csv(ab, sep='\t')
            df, inds = get_tpm_table(df,
                                     how='gene',
                                     gene_subset=gs_label,
                                     min_tpm=min_tpm,
                                     species=species)
            subset['gid'] = inds

        fname = 'figures/simplex_{}_{}.pdf'.format(source, gs_label)
        temp = ca.plot_simplex(top='splicing_ratio',
            subset=subset,
            density=True,
            density_scale=100,
            density_cmap='Purples',
            log_density=True,
            sectors=True,
            scatter=False,
            legend=False,
            title=title,
            size_scale=0.2,
            density_cbar=False,
            fname=fname)

        # add sectors from this subset to plot_df
        temp = assign_gisx_sector(temp)
        plot_df = pd.concat([temp, plot_df])

    # create the bar plots
    fname = 'figures/{}_genes_per_sector.pdf'.format(gs_label)
    make_sector_source_bar_plots(plot_df, fname, species)
    
def plot_major_iso_simplex(h5, gene):
    ca = cerberus.read(h5)
    
    # plotting settings
    c_dict, order = get_biosample_colors()
    c_dict[np.nan] = 'k'
    mm_tissues = ['adrenal gland', 'heart',
                  'muscle', 'brain', 'pgp1_excite_neuron',
                  'pgp1_astro', 'h9_osteocyte']
    mmap = {'v40': '*', 'v29': 'x', 'obs_det': '^', 'cerberus': '2', 'obs_major': '^'}
    subset = {'source': ['v40', 'obs_det', 'sample_major']}
    
    fname = f'figures/simplex_{gene.lower()}_major.pdf'
    
    df = ca.plot_simplex(gene=gene,
            hue='sample',
            cmap=c_dict,
            size='gene_tpm', 
            log_size=True,
            sectors=True,
            marker_style='source',
            mmap=mmap,
            legend=False,
            jitter=True,
            subset={'source': ['v40', 'obs_det', 'sample_major']},
            size_scale=0.2,
            fname=fname)    
    
def plot_obs_obs_major_density_simplex(h5, gene):
    
    ca = cerberus.read(h5)
    
    # sample_det
    c_dict, order = get_biosample_colors()
    c_dict[np.nan] = 'k'
    mmap = {'v40': '*', 'v29': 'x', 'obs_det': '^', 'cerberus': '2', 'sample_det_centroid': 's'} 
    subset = {'source': ['v40', 'obs_det', 'sample_det']}
    fname=f'figures/simplex_{gene.lower()}_det.pdf'
    df = ca.plot_simplex(top='splicing_ratio', 
                gene=gene,
                hue='sample',
                cmap=c_dict,
                density=True,
                density_scale=50,
                density_cmap='Purples',
                size='gene_tpm', 
                log_size=True,
                sectors=True,
                marker_style='source',
                mmap=mmap,
                legend=False,
                jitter=True,
                subset=subset,
                size_scale=0.2,
                fname=fname)
    
    
    # sample_major
    mmap = {'v40': '*', 'v29': 'x', 'obs_major': '^', 'cerberus': '2', 'sample_det_centroid': 's'} 
    subset = {'source': ['v40', 'obs_major', 'sample_major']}
    fname=f'figures/simplex_{gene.lower()}_major.pdf'
    df = ca.plot_simplex(top='splicing_ratio', 
                gene=gene,
                hue='sample',
                cmap=c_dict,
                density=True,
                density_scale=50,
                density_cmap='Purples',
                size='gene_tpm', 
                log_size=True,
                sectors=True,
                marker_style='source',
                mmap=mmap,
                legend=False,
                jitter=True,
                subset=subset,
                size_scale=0.2,
                fname=fname)
    
def sector_sankey(h5,
                  source1, source2,
                  ofile,
                  ver=None,
                  gene_subset=None):
    
    ca = cerberus.read(h5)
    
    # limit to sources
    df1 = ca.triplets.loc[ca.triplets.source==source1].copy(deep=True)
    df2 = ca.triplets.loc[ca.triplets.source==source2].copy(deep=True)
                
    # merge on gitd
    df = df1.merge(df2, how='outer', on='gid', suffixes=(f'_{source1}', f'_{source2}'))
    
    # limit to gene subset
    if gene_subset:
        gene_df, _, _ = get_gtf_info(how='gene', ver=ver)
        gene_df['gid_stable'] = cerberus.get_stable_gid(gene_df, 'gid')
        df = df.merge(gene_df[['gid_stable', 'biotype_category']],
                          how='left', left_on='gid', right_on='gid_stable')
        df = df.loc[df.biotype_category == gene_subset]
        df.drop(['biotype_category', 'gid_stable'], axis=1, inplace=True)
        
    # count numer of things
    gb_cols = [f'sector_{source1}', f'sector_{source2}']
    keep_cols = gb_cols + ['gid']
    df = df[keep_cols].groupby(gb_cols).count().reset_index()
    
    fig = plot_sankey(df,
                      source='sector_obs_det',
                      sink='sector_obs_major',
                      counts='gid',
                      color='sector',
                      title='')
    h = 700
    w = 1.8792590838529746*h
    pio.write_image(fig, ofile, width=w, height=h)
    
    return df

def plot_sr_brain_tissue_cell_line_umap(ab,
                                        lib_metadata,
                                        min_tpm,
                                        gene_subset,
                                        sample,
                                        fname):
    df = pd.read_csv(ab, sep='\t')
    df, gids = get_tpm_table(df,
                     how='sr',
                     min_tpm=min_tpm,
                     sample=sample,
                     gene_subset=gene_subset)
    df = df.transpose()

    X = df.values

    # obs table w/ sample info
    obs = df.index.to_frame()
    obs.rename({0: 'dataset'}, axis=1, inplace=True)
    obs['biosample'] = obs.dataset.str.rsplit('_', n=2, expand=True)[0]
    tissue_df = get_tissue_metadata()
    tissue_df = tissue_df[['tissue', 'biosample']]
    obs = obs.merge(tissue_df, how='left', on='biosample')
    obs.loc[obs.tissue.isnull(), 'tissue'] = obs.loc[obs.tissue.isnull(), 'biosample']
    obs.drop('biosample', axis=1, inplace=True)
    obs.rename({'tissue': 'sample'}, axis=1, inplace=True)
    obs['sample'].unique()
    obs['dataset'] = obs['dataset'].astype('string')

    # var
    var = pd.DataFrame(data=df.columns.tolist(), columns=['gid'])
    var['gid'] = var['gid'].astype('string')

    adata = anndata.AnnData(obs=obs, var=var, X=X)

    # limit only to samples that are in the long-read
    c_dict, order = get_biosample_colors()
    adata = adata[adata.obs['sample'].isin(order)]
    
    meta = pd.read_csv(lib_metadata, sep='\t')
    meta = meta[['sample', 'tissue_or_cell_line']].drop_duplicates()

    # add colors
    # map values in order specific to
    obs_col = 'sample'
    cmap, order = get_biosample_colors()
    adata.obs[obs_col] = adata.obs[obs_col].astype('category')
    obs_order = list(adata.obs_names)
    sample_order = adata.obs[obs_col].cat.categories.tolist()
    sample_colors = [cmap[s] for s in sample_order]
    adata.uns['{}_colors'.format(obs_col)] = sample_colors

    # also store rgb values in dict for use with gen_report
    for key, item in cmap.items():
        item = item[1:]
        r,g,b = tuple(int(item[i:i+2], 16) for i in (0, 2, 4))
        cmap[key] = (r,g,b)
    adata.uns['{}_dict'.format(obs_col)] = cmap

    # add tissue / cell line info
    adata.obs = adata.obs.merge(meta, how='left', on='sample')
    obs_col = 'brain_tissue_cell_line'
    adata.obs[obs_col] = adata.obs.tissue_or_cell_line
    adata.obs.loc[adata.obs['sample']=='brain', obs_col] = 'brain'

    # add colors
    # map values in order specific to
    cmap, order = get_tissue_cell_line_colors()
    cmap2, _ = get_biosample_colors()
    brain_color = cmap2['brain'] 
    order += ['brain']
    cmap['brain'] = brain_color
    adata.obs[obs_col] = adata.obs[obs_col].astype('category')
    obs_order = list(adata.obs_names)
    sample_order = adata.obs[obs_col].cat.categories.tolist()
    sample_colors = [cmap[s] for s in sample_order]
    adata.uns['{}_colors'.format(obs_col)] = sample_colors
    
    adata.obs.set_index('dataset', inplace=True)
    adata.var.set_index('gid', inplace=True)
    
    # # normalize data matrix to 10,000 count
    # sc.pp.normalize_total(adata, target_sum=1e4)

    # log 
    sc.pp.log1p(adata)
    adata.raw = adata

    # find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes = 20000, flavor = 'seurat',
        min_mean=0.0125, max_mean=3, min_disp=0.5)
    # sc.pp.highly_variable_genes(adata, n_top_genes = 5000, flavor = 'seurat',
    #     min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pp.scale(adata, max_value=10)

    # sc.tl.pca(adata, use_highly_variable=False)
    sc.tl.pca(adata, use_highly_variable=True)

    sc.pp.neighbors(adata, metric='cosine')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    
    sc.pl.umap(adata, color='brain_tissue_cell_line', frameon=True, size=120, show=False)

    f = plt.gcf()
    f.savefig(fname, dpi=500, bbox_inches='tight')

def human_v_mouse_sectors(h_h5, m_h5, h_source, m_source, 
                          h_subset=None,
                          m_subset=None,
                          gene_subset=None, 
                          h_ver='v40_cerberus',
                          m_ver='vM25_cerberus',
                          ofile='figures/human_mouse_triplets_sankey.pdf'):
    """
    Compare and plot sets of triplets between human and mouse
    
    Parameters:
        h_h5 (str): Path to human h5 file
        m_h5 (str): Path to mouse h5 file
        h_source (str): Name of source to use for human
        m_source (str): Name of source to use for mouse
        h_subset (dict): Names of samples to use in human ie {'sample': 'c2c12'}
        m_subset (dict): Names of samples to use in mouse ie {'tissue': 'cortex'}
        gene_subset (list of str): Gene biotypes to include
        h_ver (str): Name of annotation from human to pull biotypes from
        m_ver (str): Name of annotation from mouse to pull biotypes from
        
    """
    # read in triplets and limit to source
    ca = cerberus.read(h_h5)
    h_df = ca.triplets.loc[ca.triplets.source == h_source]

    ca = cerberus.read(m_h5)
    m_df = ca.triplets.loc[ca.triplets.source == m_source]    
    
    def subset_on_dict(df, subset):
        for col, val in subset.items():
                if type(val) != list:
                    val = [val]
                df = df.loc[df[col].isin(val)]
        return df
                
    # limit to samples if given
    if h_subset:
        h_df = subset_on_dict(h_df, h_subset)
    if m_subset:
        m_df = subset_on_dict(m_df, m_subset)

    # get sectors for each gene
    h_df = assign_gisx_sector(h_df)
    m_df = assign_gisx_sector(m_df)
                    
    # merge in gene info and subset
    if gene_subset:
        gene_df, _, _ = get_gtf_info(how='gene', ver=h_ver)
        gene_df['gid'] = cerberus.get_stable_gid(gene_df, 'gid')
        h_df = h_df.merge(gene_df[['gid', 'biotype_category']],
                          how='left', on='gid')
        
        gene_df, _, _ = get_gtf_info(how='gene', ver=m_ver)
        gene_df['gid'] = cerberus.get_stable_gid(gene_df, 'gid')
        m_df = m_df.merge(gene_df[['gid', 'biotype_category']],
                          how='left', on='gid')
        
        h_df = h_df.loc[h_df.biotype_category.isin(gene_subset)]
        m_df = m_df.loc[m_df.biotype_category.isin(gene_subset)]
    
    # get matching gids from human and mouse
    df = get_human_mouse_gid_table('../../refs/biomart_human_to_mouse.tsv')
    
    df = df.merge(h_df, how='right', left_on='Gene stable ID',
                  right_on='gid')
        
    df = df.merge(m_df, how='right', left_on='Mouse gene stable ID',
                  right_on='gid', suffixes=('_human', '_mouse'))
    
    # get only genes that are both expressed
    df = df.loc[~(df.sector_human.isnull())&~(df.sector_mouse.isnull())]
    
    # count up instances of each sector pair
    df = df[['Gene stable ID', 'sector_human', 'sector_mouse']].groupby(['sector_human', 'sector_mouse']).count().reset_index()
    df.rename({'Gene stable ID': 'n_genes'}, axis=1, inplace=True)                  
    
    # display some stuff
    print('Total number of human genes / sector')
    print(df.groupby(['sector_human']).sum().reset_index()[['sector_human', 'n_genes']])

    print('Total number of mouse genes / sector')
    print(df.groupby(['sector_mouse']).sum().reset_index()[['sector_mouse', 'n_genes']])
    
    print('Total number of conserved genes / sector')
    temp = df.loc[df.sector_human==df.sector_mouse].groupby('sector_human').sum().reset_index()[['sector_human', 'n_genes']]
    temp2 = df.groupby(['sector_human']).sum().reset_index()[['sector_human', 'n_genes']]
    temp2.rename({'n_genes': 'n_total_genes'}, axis=1, inplace=True)
    temp = temp.merge(temp2, on='sector_human')
    temp['perc'] = (temp.n_genes/temp.n_total_genes)*100
    print(temp)

    fig = plot_sankey(df,
                  source='sector_human',
                  sink='sector_mouse',
                  counts='n_genes',
                  color='sector',
                  title='')
    pio.write_image(fig, ofile, width=700, height=700)
    
    n_genes = df.n_genes.sum(axis=0)
    n_cons_genes = df.loc[df.sector_human == df.sector_mouse].n_genes.sum(axis=0)
    print('{}/{} ({:.2f}%) of genes have conserved sectors between human and mouse.'.format(n_cons_genes, n_genes, (n_cons_genes/n_genes)*100))
    
    return df

def plot_human_mouse_simplex(ca, m_ca, h_gene, m_gene, odir):
    

    # human
    c_dict, order = get_biosample_colors()
    c_dict[np.nan] = 'k'
    mmap = {'v40': '*', 'v29': 'x', 'obs_mm_det': '^', 'cerberus': '2', 'sample_mm_det_centroid': 'x'} 
    subset = {'source': ['v40', 'obs_mm_det', 'sample_mm_det', 'sample_mm_det_centroid']}
    
    df = ca.plot_simplex(top='splicing_ratio', 
                gene=h_gene,
                hue='sample',
                cmap=c_dict,
                size='gene_tpm', 
                log_size=True,
                sectors=True,
                marker_style='source',
                mmap=mmap,
                legend=False,
                jitter=True,
                subset=subset,
                size_scale=0.5,
                fname='{}/simplex_{}.pdf'.format(odir, h_gene.lower()))

    # mouse
    c_dict, order = get_lr_bulk_sample_colors()
    c_dict[np.nan] = 'k'
    mmap = {'vM25': '*', 'v29': 'x', 'obs_det': '^', 'cerberus': '2', 'sample_det_mouse_centroid': 'x'} 
    subset = {'source': ['vM25', 'obs_det', 'sample_det_mouse', 'sample_det_mouse_centroid']}
    df = m_ca.plot_simplex(top='splicing_ratio', 
                gene=m_gene,
                hue='sample',
                cmap=c_dict,
                size='gene_tpm', 
                log_size=True,
                sectors=True,
                marker_style='source',
                mmap=mmap,
                legend=False,
                jitter=True,
                subset=subset,
                size_scale=0.5,
                fname='{}/simplex_mouse_{}.pdf'.format(odir, m_gene.lower()))
    
def plot_n_transcripts_per_sample(ca,
                                  gene,
                                  source,
                                  fname):

    temp = ca.triplets.loc[(ca.triplets.source==source)&\
                           (ca.triplets.gname==gene)].sort_values(by='n_iso', ascending=False)
    
    c_dict, order = get_biosample_colors()
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    plt.figure(figsize=(4,4))

    order = temp['sample'].tolist()
    ax = sns.barplot(temp, y='n_iso', x='sample', palette=c_dict, order=order, saturation=1)
    _ = plt.xticks([], [])
    _ = ax.set(xlabel='Sample', ylabel=f'# {gene} transcripts')

    # add hline at 1
    xlims = ax.get_xlim()
    color = '#5c5c5c'
    ax.hlines(1, xlims[0], xlims[1],
                  colors=color, linestyles='dashed',
                  linewidth=2)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.savefig(fname, dpi=500, bbox_inches='tight')
    
def get_counts_per_gene(filt_ab, obs_col, feat, min_tpm, gene_subset):
    # counts per gene
    df = pd.read_csv(filt_ab, sep='\t')
    df = get_det_table(df,
                        groupby=obs_col,
                        how=feat, 
                        gene_subset=gene_subset,
                        min_tpm=min_tpm)

    df = df.transpose()
    df.reset_index(inplace=True)
    df.rename({'index':'feat_id'}, axis=1, inplace=True)
    if feat in ['tss', 'ic', 'tes']:
        df['gid'] = df.feat_id.str.split('_', expand=True)[0]
    else:
        df['gid'] = df.feat_id.str.split('[', expand=True)[0]
    df = df[['gid', 'feat_id']].groupby('gid').count().reset_index()
    df.rename({'feat_id': 'counts'}, axis=1, inplace=True)
    
    return df

def get_counts_per_sample_per_gene(filt_ab, obs_col, feat, min_tpm, gene_subset):
    
    # counts per sample per gene
    df = pd.read_csv(filt_ab, sep='\t')
    df = get_det_table(df,
                        groupby=obs_col,
                        how=feat, 
                        gene_subset=gene_subset,
                        min_tpm=min_tpm)

    df = df.melt(var_name='feat_id', value_name='det', ignore_index=False)
    df.reset_index(inplace=True)
    df.rename({'index': obs_col}, axis=1, inplace=True)
    df = df.loc[df.det == True]
    df.drop('det', axis=1, inplace=True)

    if feat in ['tss', 'ic', 'tes']:
        df['gid'] = df.feat_id.str.split('_', expand=True)[0]
    else:
        df['gid'] = df.feat_id.str.split('[', expand=True)[0]
    df = df.groupby(['biosample', 'gid']).count().reset_index()
    df.rename({'feat_id': 'counts'}, axis=1, inplace=True)
    
    return df

def plot_feats_per_sample_gene(temp, feat, obs_col):
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42

    if feat in ['tss', 'tes', 'ic']:
        x_hr = feat.upper()
    elif feat == 'iso':
        x_hr = 'isoform'
    c_dict, order = get_feat_triplet_colors_2(feat)
    color = c_dict[feat]
    
    # get each df
    temp1 = temp.loc[temp.source == 'Sample']
    temp1, order, col = renum_max_feats(temp1, 'counts', 10)
    temp1 = temp1[[col, 'biosample', 'gid']].groupby([col, 'biosample']).count().reset_index()
    temp1.rename({'gid': 'n_feats'}, axis=1, inplace=True)
    temp1_max = temp1['n_feats'].max()
    
    temp2 = temp.loc[temp.source == 'Observed']
    temp2, order, col = renum_max_feats(temp2, 'counts', 10)
    temp2 = temp2[[col, 'gid']].groupby([col]).count().reset_index()
    temp2.rename({'gid': 'n_feats'}, axis=1, inplace=True)
    temp2_max = temp2['n_feats'].max()
    
    ymax = max(temp1_max, temp2_max)
    ylim = (0, ymax)
    
    # feats / sample / gene
    ax = sns.catplot(data=temp1, y='n_feats', x=col, kind='bar',
                     edgecolor=None, saturation=1,
                     alpha=1, color=color, order=order)
    xlabel = f'# {x_hr}s / gene / {obs_col}'
    ylabel = '# genes'
    ax.set(xlabel=xlabel, ylabel=ylabel, ylim=ylim)
    fname = f'{fig_dir}/{feat}_gene_sample_hist.pdf'
    plt.savefig(fname, dpi=500)

    # observed overall feats / gene
    plt.figure(figsize=(3,4))
    ax = sns.catplot(data=temp2, y='n_feats', x=col, kind='bar',
                     edgecolor=None, saturation=1,
                     alpha=1, color=color, order=order)
    xlabel = f'# {x_hr}s / gene'
    ylabel = '# genes'
    ax.set(xlabel=xlabel, ylabel=ylabel, ylim=ylim)
    fname = f'{fig_dir}/{feat}_gene_hist.pdf'
    plt.savefig(fname, dpi=500)

def plot_feats_per_sample_gene_ecdf(temp, feat, obs_col):
    from matplotlib.lines import Line2D

    c_dict, order = get_feat_triplet_colors_2(feat)
    color = c_dict[feat]

    plt.figure(figsize=(4,4))
    temp1 = temp.loc[temp.source=='Observed']
    ax = sns.ecdfplot(data=temp1, x='counts', color=color, linewidth=3, linestyle='--')
    temp1 = temp.loc[temp.source=='Sample']
    ax = sns.ecdfplot(data=temp1, x='counts', color=color, linewidth=3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    legend_elements = [Line2D([0], [0], color=color, lw=3, label='Sample'),
                       Line2D([0], [0], color=color, lw=3, linestyle='--', label='Observed')]
    ax.legend(handles=legend_elements)

    x_hr = feat.upper()
    if feat in ['tss', 'ic', 'tes']:
        xlabel_thing = f'# {x_hr}s'
    elif feat == 'iso':
        xlabel_thing = '# transcripts'
    xlabel = f'{xlabel_thing} / gene'
    _ = ax.set(xlabel=xlabel)

    fname = f'{fig_dir}/{feat}_gene_sample_ecdf.pdf'
    plt.savefig(fname, dpi=500)
    
def renum_max_feats(temp, col, max_val):
    """
    Parameters:
        temp (df): DataFrame
        col (str): Col to transform
        max_val (int): Max value
    """
    max_str = '{}+'.format(max_val)
    temp['rank'] = temp[col].astype(int).astype(str)
    temp.loc[temp[col]>=max_val, 'rank'] = max_str
    order = [str(int(i)) for i in range(1,max_val)]
    order += [max_str]
    return temp, order, 'rank'

def perc(n_num, n):
    return (n_num/n)*100

def plot_perc_one_feat(df, feat, opref='figures/'):
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['pdf.fonttype'] = 42
    plt.figure(figsize=(2,4))

    c_dict, order = get_feat_triplet_colors_2(feat)
    color = c_dict[feat]
    
    if feat in ['tss', 'ic', 'tes']:
        ylabel = f'% genes with 1 {feat.upper()}'
    elif feat == 'iso':
        ylabel = '% genes with 1 isoform'
    order = ['Gene', 'Gene / sample']
    ax = sns.barplot(df, y='perc', x='combo', color=color, order=order, saturation=1)
    ax.set(ylim=(0, 100), ylabel=ylabel, xlabel='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    _ = plt.setp(ax.get_xticklabels(), rotation=90)
    
    fname = f'{opref}/perc_genes_one_{feat}.pdf'
    plt.savefig(fname, dpi=500, bbox_inches='tight')
    
def get_n_samples_vs_n_feats(df):
    """
    Compute how many samples each feature is expressed in
    """
    temp = df.copy(deep=True)
    
    # compute how many samples each gene is expressed in 
    temp1 = temp.loc[temp.source=='Sample']
    temp1 = temp1[['gid', 'biosample']].groupby('gid').count().reset_index()
    temp1.rename({'biosample':'n_samples'}, axis=1, inplace=True)
    temp1.head()
    
    temp2 = temp.loc[temp.source=='Observed'].drop('biosample', axis=1)
    temp2 = temp2.merge(temp1, on='gid', how='left')
    temp2.head()
    
    return temp2

def plot_n_samples_vs_n_feats(temp, feat, opref):
    c_dict, order = get_feat_triplet_colors_2(feat)
    color = c_dict[feat]
    
    df = get_n_samples_vs_n_feats(temp)
    # ax = sns.jointplot(df, x='n_samples', y='counts', color=color, marginal_kws={'discrete':True, 'edgecolor':None}, joint_kws={'alpha':0.25, 's':10})
    ax = sns.jointplot(df, x='n_samples', y='counts', kind='hex', color=color, marginal_kws={'discrete':True, 'edgecolor':None},
                       joint_kws={'edgecolors':'none', 'bins':'log', 'mincnt':0})

    ylabel = '# observed '
    if feat in ['tss', 'ic', 'tes']:
        ylabel += f'{feat.upper()}s'
    elif feat == 'iso':
        ylabel += 'transcripts'
    xlabel = '# samples'
    
    ax.ax_joint.set(xlabel=xlabel, ylabel=ylabel)
    
    fname = f'{opref}/n_{feat}s_per_n_samples_expressed.pdf'
    plt.savefig(fname, dpi=500, bbox_inches='tight')
    
def plot_n_predom_transcripts(pi_tpm_file,
                              ver,
                              gene_subset,
                              fname,
                              max_isos=None,
                              figsize=(6,6)):
    df = pd.read_csv(pi_tpm_file, sep='\t')
    
    # only predominant transcripts
    df = df.loc[df.triplet_rank==1]

    # count number of unique predominant transcripts
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

    return df