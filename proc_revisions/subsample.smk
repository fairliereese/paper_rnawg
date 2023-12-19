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
    temp = read_annot_to_counts(df)
    temp = add_ab_metdata(temp, talon_ab)

    temp.to_csv(ofile, sep='\t', index=False)

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

sample_depths = [.1, .15, .2, .4, .6, .8, .9, 1]

rule subsample_wtc11_transcript:
    input:
        ab = config['lr']['cerberus']['filt_ab']
    resources:
        mem_gb = 32,
        threads = 2
    output:
        ab = config['lr']['subsample']['filtab']
    run:
        subsample_ab(input.ab,
                     wildcards.depth,
                     output.ab)

rule subsample_wtc11_transcript_summary:
    input:
        files = expand(config['lr']['subsample']['filt_ab'],
               species='human',
               subsample_depth=subsample_depths)
    resources:
        mem_gb = 32,
        threads = 2
    params:
        min_tpm = 1,
        gene_subset = 'polya'
    output:
        ofile = config['lr']['subsample']['transcript_summary']
    run:
        df = pd.DataFrame()
        df['file'] = input.files
        df['depth'] = df.file.str.rsplit('_', n=1, expand=True)[0]

        n_genes = []
        for ind, entry in df.iterrows():
            file = entry.file
            temp = pd.read_csv(file, sep='\t')
            _, inds = get_tpm_table(temp,
                    how='iso',
                    gene_subset=params.gene_subset,
                    min_tpm=params.min_tpm)
            n_genes += len(inds)
        df['n_transcripts'] = n_genes
        df.to_csv(ofile, sep='\t', index=False)

rule subsample_wtc11_gene:
    input:
        read_annot = config['lr']['talon']['full_annot'],
        ab = config['lr']['talon']['fusion_fix']['ab']
    resources:
        mem_gb = 32,
        threads = 2
    output:
        ab = config['lr']['subsample']['ab']
    run:
        subsample_read_annot(input.read_annot,
                             input.ab,
                             wildcards.depth,
                             output.ab)

rule subsample_wtc11_gene_summary:
    input:
        files = expand(config['lr']['subsample']['ab'],
               species='human',
               subsample_depth=subsample_depths)
    resources:
        mem_gb = 32,
        threads = 2
    params:
        min_tpm = 1,
        gene_subset = 'polya'
    output:
        ofile = config['lr']['subsample']['gene_summary']
    run:
        df = pd.DataFrame()
        df['file'] = input.files
        df['depth'] = df.file.str.rsplit('_', n=1, expand=True)[0]

        n_genes = []
        for ind, entry in df.iterrows():
            file = entry.file
            temp = pd.read_csv(file, sep='\t')
            _, inds = get_tpm_table(temp,
                    how='gene',
                    gene_subset=params.gene_subset,
                    min_tpm=params.min_tpm)
            n_genes += len(inds)
        df['n_genes'] = n_genes
        df.to_csv(ofile, sep='\t', index=False)

rule all_subsample:
    input:
        expand(config['lr']['subsample']['gene_summary'],
               species='human'),
        expand(config['lr']['subsample']['transcript_summary'],
               species='human')
