sample_depths = [.1, .15, .2, .4, .6, .8, .9]
sample_reps = [i for i in range(10)]
sample_gene_subsets = ['polya', 'protein_coding', 'pseudogene', 'lncRNA']

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

    depth = float(depth)

    # read in the transcript id, dataset, and the read name from the read annot
    df = pd.read_csv(read_annot, sep='\t',
                     usecols=[0,1,12])

    # limit to just the wtc11 datasets
    wtc11_cols = ['wtc11_1_3', 'wtc11_1_2', 'wtc11_1_1']
    df = df.loc[df.dataset.isin(wtc11_cols)]

    # sample to target depth
    sample_df = df.sample(frac=depth,
                          replace=False)
    temp = read_annot_to_counts(sample_df)
    temp = add_ab_metdata(temp, talon_ab)

    temp.to_csv(ofile, sep='\t', index=False)

def subsample_ab(talon_ab,
                 depth,
                 ofile):
    """
    Subsample WTC11 data to a specified depth from the abundance matrix
    and output a TALON abundance file.
    """

    depth = float(depth)

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

rule subsample_wtc11_transcript:
    input:
        ab = config['lr']['cerberus']['filt_ab']
    resources:
        mem_gb = 32,
        threads = 2
    output:
        ab = config['lr']['subsample']['filt_ab']
    run:
        subsample_ab(input.ab,
                     wildcards.subsample_depth,
                     output.ab)

rule subsample_wtc11_transcript_summary:
    input:
        files = expand(config['lr']['subsample']['filt_ab'],
               species='human',
               subsample_depth=sample_depths,
               subsample_rep=sample_reps)
    resources:
        mem_gb = 32,
        threads = 2
    params:
        min_tpm = 0
    output:
        ofile = config['lr']['subsample']['transcript_summary']
    run:
        df = pd.DataFrame()
        df['file'] = list(input.files)
        df['depth'] = df.file.str.rsplit('_', n=2, expand=True)[1]
        df['rep'] = df.file.str.rsplit('_', n=1, expand=True)[1].str.rsplit('.', n=1, expand=True)[0]

        files_loop = df.file.tolist()
        depths_loop = df.depth.tolist()
        reps_loop = df.rep.tolist()

        n_genes = []
        files = []
        depths = []
        reps = []
        subsets = []

        for file, depth, rep in zip(files_loop, depths_loop, reps_loop):
            temp = pd.read_csv(file, sep='\t')
            for gene_subset in sample_gene_subsets:
                _, inds = get_tpm_table(temp,
                        how='iso',
                        gene_subset=gene_subset,
                        min_tpm=params.min_tpm)
                n_genes.append(len(inds))
                files.append(file)
                depths.append(depth)
                reps.append(rep)
                subsets.append(gene_subset)

        df = pd.DataFrame()
        df['n_transcripts'] = n_genes
        df['depth'] = depths
        df['rep'] = reps
        df['gene_subset'] = subsets

        df.to_csv(output.ofile, sep='\t', index=False)

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
                             wildcards.subsample_depth,
                             output.ab)

def get_corr_t_summary(files, full_ab, ofile, params):

    file_df = pd.DataFrame()
    file_df['file'] = list(files)
    file_df['depth'] = file_df.file.str.rsplit('_', n=2, expand=True)[1]
    file_df['rep'] = file_df.file.str.rsplit('_', n=1, expand=True)[1].str.rsplit('.', n=1, expand=True)[0]

    spearman_corrs = []
    pearson_corrs = []
    gs = []
    reps = []
    depths = []

    gene_subsets = params['gene_subsets']
    min_tpm = params['min_tpm']

    for ind, entry in file_df.iterrows():
        file = entry.file
        depth = entry.depth
        rep =  entry.depth
        for g in gene_subsets:
            ab_df = pd.read_csv(full_ab, sep='\t')
            ab_df, _ = get_tpm_table(ab_df,
                                 how='iso',
                                 gene_subset=g,
                                 min_tpm=min_tpm,
                                  groupby='sample',
                                 sample=['wtc11'])
            ab_df = ab_df.reset_index()
            ab_df.rename({'index': 'tid'}, axis=1, inplace=True)

            sub_ab_df = pd.read_csv(file, sep='\t')
            sub_ab_df, _ = get_tpm_table(sub_ab_df,
                                 how='iso',
                                 gene_subset=g,
                                 min_tpm=min_tpm,
                                      groupby='sample',
                                 sample=['wtc11'])
            sub_ab_df = sub_ab_df.reset_index()
            sub_ab_df.rename({'index': 'tid'}, axis=1, inplace=True)

            # intersection (only correlate stuff detected in both)
            temp = ab_df.merge(sub_ab_df, how='inner',
                  on='tid',
                  suffixes=('_full', '_subs'))

            x = 'wtc11_full'
            y = 'wtc11_subs'
            try:
                rho, p = st.spearmanr(temp[x].tolist(), temp[y].tolist())
                r, p2 = st.pearsonr(temp[x].tolist(), temp[y].tolist())
            except:
                import pdb; pdb.set_trace()
            spearman_corrs.append(r)
            pearson_corrs.append(rho)
            gs.append(g)
            reps.append(rep)
            depths.append(depth)

    corr_df = pd.DataFrame()
    corr_df['depth'] = depths
    corr_df['rep'] = reps
    corr_df['gene_subset'] = gs
    corr_df['pearson_r'] = pearson_corrs
    corr_df['spearman_rho'] = spearman_corrs
    corr_df.to_csv(ofile, sep='\t', index=False)

rule subsample_wtc11_transcript_corr_summary:
    input:
        files = expand(config['lr']['subsample']['filt_ab'],
               species='human',
               subsample_depth=sample_depths,
               subsample_rep=sample_reps),
        filt_ab = expand(config['lr']['cerberus']['filt_ab'],
                species='human')[0]
    resources:
        mem_gb = 32,
        threads = 2
    params:
        min_tpm = 0,
        gene_subsets = ['protein_conding', 'lncRNA', 'pseudogene']
    output:
        ofile = config['lr']['subsample']['transcript_corr_summary']
    run:
        get_corr_t_summary(input.files, input.filt_ab, output.ofile, params)

rule subsample_wtc11_gene_summary:
    input:
        files = expand(config['lr']['subsample']['ab'],
               species='human',
               subsample_depth=sample_depths,
               subsample_rep=sample_reps)
    resources:
        mem_gb = 32,
        threads = 2
    params:
        min_tpm = 0
    output:
        ofile = config['lr']['subsample']['gene_summary']
    run:
        df = pd.DataFrame()
        df['file'] = list(input.files)
        df['depth'] = df.file.str.rsplit('_', n=2, expand=True)[1]
        df['rep'] = df.file.str.rsplit('_', n=1, expand=True)[1].str.rsplit('.', n=1, expand=True)[0]

        files_loop = df.file.tolist()
        depths_loop = df.depth.tolist()
        reps_loop = df.rep.tolist()

        n_genes = []
        files = []
        depths = []
        reps = []
        subsets = []

        for file, depth, rep in zip(files_loop, depths_loop, reps_loop):
            temp = pd.read_csv(file, sep='\t')
            for gene_subset in sample_gene_subsets:
                _, inds = get_tpm_table(temp,
                            how='gene',
                            gene_subset=gene_subset,
                            min_tpm=params.min_tpm)
                n_genes.append(len(inds))
                files.append(file)
                depths.append(depth)
                reps.append(rep)
                subsets.append(gene_subset)

        df = pd.DataFrame()
        df['n_genes'] = n_genes
        df['depth'] = depths
        df['rep'] = reps
        df['gene_subset'] = subsets

        df.to_csv(output.ofile, sep='\t', index=False)

# rule subsample_calc_triplets:
#     input:
#         h5 = config['lr']['cerberus']['ca_annot'],
#         filt_ab = config['lr']['subsample']['filt_ab']
#     params:
#         min_tpm = 0,
#         gene_subset = 'polya'
#     resources:
#         mem_gb = 16,
#         threads = 1
#     output:
#         h5 = config['lr']['subsample']['ca_triplets']
#     run:
#         ca = cerberus.read(input.h5)
#         filt_ab_df = pd.read_csv(input.filt_ab, sep='\t')
#
#         # observed triplets
#         source = f'{wildcards.subsample_depth}_obs_det'
#         df, tids = get_tpm_table(filt_ab_df,
#                    how='iso',
#                    min_tpm=params.min_tpm,
#                    subset=params.gene_subset)
#         df = ca.get_subset_triplets(tids, source)
#         ca.add_triplets(df)
#         ca.write(output.h5)

rule all_subsample:
    input:
        expand(config['lr']['subsample']['gene_summary'],
               species='human'),
        expand(config['lr']['subsample']['transcript_summary'],
               species='human'),
        expand(config['lr']['subsample']['transcript_corr_summary'],
               species='human')
        # expand(config['lr']['subsample']['ca_triplets'],
        #        species='human',
        #        subsample_depth=sample_depths,
        #        subsample_rep=sample_reps)
