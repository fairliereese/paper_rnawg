def fix_talon_fusion_transcripts(talon_ab,
                                 talon_filt_ab,
                                 talon_gtf,
                                 ref_ics,
                                 ref_gtf,
                                 wc,
                                 ofile_gtf,
                                 ofile_ab):
    """
    Fix fusion gene assignments based on the ss concordance between
    annotated (v29 or vM21) transcripts. Update gene counts accordingly.
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
    
    # df, gene_id == new gene id, Name == transcript name
    # also update the unfiltered abundance table
    temp = df[['Name', 'gene_id']]
    temp.rename({'Name': 'annot_transcript_id',
                 'gene_id_known': 'gid_stable'},
                axis=1, inplace=True)
    df = pd.read_csv(talon_ab, sep='\t')
    df['gid_stable'] = cerberus.get_stable_gid(df, 'annot_gene_id')
    gene_cols = ['gene_ID', 'annot_gene_id', 'annot_gene_name', 'gene_novelty', 'gid_stable']
    gene_info = df[gene_cols].drop_duplicates()
    temp = temp.merge(gene_info, how='left', on='gid_stable')
    temp.drop('gid_stable', axis=1, inplace=True)    
    df = df.merge(temp, how='left', on='annot_transcript_id', suffixes=('', '_fix'))
    fix_inds = df.loc[df.gene_ID_fix.notnull()].index
    gene_cols = [g for g in gene_cols if g!='gid_stable']
    for g in gene_cols:
        fix_col = f'{g}_fix'
        df.loc[fix_inds, g] = df.loc[fix_inds, fix_col]
        df.drop(fix_col, axis=1, inplace=True)
    df.to_csv(ofile_ab, sep='\t', index=False)
    
def fix_talon_known_fusion_transcripts(talon_ab,
                                       talon_filt_ab,
                                       talon_gtf,
                                       ref_ics,
                                       ref_gtf,
                                       wc,
                                       ofile_gtf,
                                       ofile_ab): 
    if wc['species'] == 'human':
        ver = 'v40_cerberus'
    elif wc['species'] == 'mouse':
        ver = 'vM25_cerberus'

    gtf_df, _, _ = get_gtf_info(ver=ver, how='gene', add_stable_gid=True)
    gtf_df = gtf_df[['gid_stable', 'gname']]
    gtf_df.head()

    talon_df = pd.read_csv(talon_filt_ab, sep='\t')
    talon_df['gid'] = cerberus.get_stable_gid(talon_df, 'annot_gene_id')

    tids = talon_df.loc[(talon_df.transcript_novelty!='Known')&(talon_df.gene_novelty=='Known'), 'annot_transcript_id'].tolist()

    df = pr.read_gtf(talon_gtf, rename_attr=True, duplicate_attr=True)
    df = cerberus.get_ic(df)
    df.rename({'transcript_id':'Name', 'ic': 'Coordinates'}, axis=1, inplace=True)
    df['source'] = 'lapa'
    
    # get gids from annotated fusion genes
    rt_df, _, _ = get_gtf_info(ver=ver, how='iso')
    rt_df = rt_df.loc[rt_df.readthrough_transcript==True]
    rt_df['gid_stable'] = cerberus.get_stable_gid(rt_df, 'gid')
    fusion_gids = rt_df.gid_stable.tolist()
    
    # get reference transcript id to ic name map :'(
    ref_df = pr.read_gtf(talon_gtf, rename_attr=True, duplicate_attr=True)
    ref_df = cerberus.get_ic(ref_df)
    ref_df.head()
    
    known_df = pd.read_csv(ref_ics, sep='\t')
    known_df['gene_id'] = known_df.Name.str.split('_', expand=True)[0]
    known_df['source'] = 'lapa'
    
    ref_df.rename({'ic':'Coordinates'}, axis=1, inplace=True)
    print(len(ref_df.index))
    ref_df= ref_df.merge(known_df[['Coordinates', 'gene_id', 'Chromosome', 'Name']], 
                 how='left', 
                 on=['Coordinates', 'gene_id', 'Chromosome'])
    print(len(ref_df.index))
    
    known_df = pd.read_csv(ref_ics, sep='\t')
    known_df['gene_id'] = known_df.Name.str.split('_', expand=True)[0]
    known_df['source'] = 'lapa'

    _, f_ss_ic_df = get_ss_from_ic(df, ['lapa'])
    f_ss_ic_df.drop(['source', 'novelty'], axis=1, inplace=True)

    _, k_ss_ic_df = get_ss_from_ic(known_df, ['lapa'])
    k_ss_ic_df.drop(['source', 'novelty'], axis=1, inplace=True)
    
    f_ss_ic_df['ss_id'] = f_ss_ic_df.Name+f_ss_ic_df.Chromosome+f_ss_ic_df.Strand+f_ss_ic_df.Start.astype(str)
    k_ss_ic_df['ss_id'] = k_ss_ic_df.Name+k_ss_ic_df.Chromosome+k_ss_ic_df.Strand+k_ss_ic_df.Start.astype(str)
    
    f_ss_ic_df = f_ss_ic_df.merge(k_ss_ic_df, how='left', on=['Chromosome', 'Strand', 'Start', 'ss_type'],
                          suffixes=('', '_known'))
    f_ss_ic_df['known'] = False
    f_ss_ic_df.loc[f_ss_ic_df.Name_known.notnull(), 'known'] = True
    
    # 1. get all transcript + gene combinations
    df = f_ss_ic_df[['gene_id', 'Name', 'gene_id_known', 'known']].drop_duplicates().copy(deep=True)
    
    # 2. compute number of splice sites / transcript
    temp = f_ss_ic_df[['Name', 'ss_id']].groupby('Name').nunique()
    temp = temp.reset_index().rename({'ss_id': 'n_ss'}, axis=1)
    df = df.merge(temp, how='left', on='Name')
    
    # 3. add in the novelty of the talon gene and transcript
    df = df.merge(talon_df[['transcript_novelty', 'gene_novelty', 'annot_transcript_id']], 
                  how='left', 
                  left_on='Name',
                  right_on='annot_transcript_id')
    df.drop('annot_transcript_id', axis=1, inplace=True)
    
    # 4. subset for novel transcripts that belong to annotated genes
    df = df.loc[(df.transcript_novelty!='Known')&\
                (df.gene_novelty=='Known')]
    
    # 5. get all transcripts that have >= 1 splice site shared w/ a readthrough gene
    tids = df.loc[df.gene_id_known.isin(fusion_gids), 'Name'].unique().tolist()
    df = df.loc[df.Name.isin(tids)]
    
    # 6. get how many unique genes share splice sites with this transcript
    temp = df.loc[df.known==True].copy(deep=True)
    temp = temp[['Name', 'gene_id_known']].groupby('Name').nunique().reset_index()
    temp.rename({'gene_id_known': 'n_genes'}, axis=1, inplace=True)
    df = df.merge(temp, how='left', on='Name')
    
    # 7. get how many unique genes share splice sites with this transcript
    temp = df.loc[df.known==True].copy(deep=True)
    temp = temp.loc[~temp.gene_id_known.isin(fusion_gids)]
    temp = temp[['Name', 'gene_id_known']].groupby('Name').nunique().reset_index()
    temp.rename({'gene_id_known': 'n_genes_no_rt'}, axis=1, inplace=True)
    df = df.merge(temp, how='left', on='Name')
    
    # 8. limit to transcripts w/ n_genes_no_rt < n_genes
    # and n_genes_no_rt == 1 (only one choice)
    # where the gene that was merged with is the non-rt one
    df = df.loc[(df.n_genes_no_rt==1)&(df.n_genes>df.n_genes_no_rt)]
    df = df.loc[~df.gene_id_known.isin(fusion_gids)]
    df.head()
    
    # 9. get the number of splice sites supported by the transcript w/ the most
    # shared splice sites for each transcript + annotated gene combination, 
    # then merge in based on transcript from obs. transcripts
    temp = f_ss_ic_df.copy(deep=True)
    temp = temp[['Name', 'ss_id', 'gene_id_known', 'Name_known']].drop_duplicates()
    temp = temp.groupby(['Name', 'gene_id_known', 'Name_known']).nunique().reset_index()
    temp.rename({'ss_id': 'n_supp_ss',
                 'Name_known': 'Name_supp'}, axis=1, inplace=True)
    temp = temp.sort_values(by='n_supp_ss', ascending=False)
    temp = temp.drop_duplicates(subset=['Name', 'gene_id_known'], keep='first')
    df = df.merge(temp, how='left', on=['Name', 'gene_id_known'])
    
    # 9.5. get the number of splice sites supported by the transcript w/ the most
    # shared splice sites for each transcript + annotated gene combination, 
    # then merge in based on transcript from obs. transcripts
    temp = f_ss_ic_df.copy(deep=True)
    temp = temp[['Name', 'ss_id', 'gene_id_known', 'Name_known']].drop_duplicates()
    temp = temp.loc[temp.gene_id_known.isin(fusion_gids)]
    temp = temp.groupby(['Name', 'gene_id_known', 'Name_known']).nunique().reset_index()
    temp.rename({'ss_id': 'n_supp_ss_rt',
                 'Name_known': 'Name_supp_rt',
                 'gene_id_known': 'gene_id_known_rt'}, axis=1, inplace=True)
    temp = temp.sort_values(by='n_supp_ss_rt', ascending=False)
    temp = temp.drop_duplicates(subset=['Name'], keep='first')
    df = df.merge(temp, how='left', on=['Name'])
    
    # 10. remove entries where talon gene id is already
    # the other choice
    df['gid'] = cerberus.get_stable_gid(df, 'gene_id')
    df = df.loc[df.gid!='gene_id_known']
    
    # 11. get the total # of splice sites per annotated transcript
    temp = k_ss_ic_df[['Name', 'ss_id']].groupby('Name').nunique().reset_index()
    temp = temp.rename({'ss_id': 'n_supp_total_ss',
                                      'Name': 'Name_supp'}, axis=1)
    df = df.merge(temp, how='left', on='Name_supp')
    temp = temp.rename({'n_supp_total_ss': 'n_supp_total_ss_rt',
                        'Name_supp': 'Name_supp_rt'}, axis=1)
    df = df.merge(temp, how='left', on='Name_supp_rt')
    
    # 12. limit to obs. transcripts that were annotated to readthrough genes
    df = df.loc[df.gid.isin(fusion_gids)]
    len(df.index)
    
    df['perc_supp_ss'] = (df.n_supp_ss/df.n_ss)*100
    df['perc_supp_ss_rt'] = (df.n_supp_ss_rt/df.n_ss)*100
    df['perc_supp_annot_ss'] = (df.n_supp_ss/df.n_supp_total_ss)*100
    df['perc_supp_annot_ss_rt'] = (df.n_supp_ss_rt/df.n_supp_total_ss_rt)*100
    
    temp = df.loc[df.gid.isin(fusion_gids)]
    len(temp.Name.unique())
    
    df['perc_supp_diff'] = df['perc_supp_ss']-df['perc_supp_ss_rt']
    df['perc_supp_annot_diff'] = df['perc_supp_annot_ss']-df['perc_supp_annot_ss_rt']
    df = df.loc[df.perc_supp_annot_ss.notnull()]
    
    # 14. compute distances from transcript start and end to closest start and end from any transcript in
    talon_df = pr.read_gtf(talon_gtf, rename_attr=True, duplicate_attr=True)
    talon_df = talon_df.df
    talon_df = talon_df.loc[talon_df.Feature=='transcript']
    
    known_df = pr.read_gtf(ref_gtf, rename_attr=True, duplicate_attr=True)
    known_df = known_df.df
    known_df = known_df.loc[known_df.Feature=='transcript']
    
    known_df['gene_id'] = cerberus.get_stable_gid(known_df, 'gene_id')
    # add coords to novel thing (df)
    temp = talon_df[['transcript_id', 'Start', 'End', 'Chromosome', 'Strand']]
    temp.rename({'transcript_id':'Name'}, axis=1, inplace=True)
    df = df.merge(temp, how='left', on='Name')
    
    df = df.drop_duplicates()
    temp3 = pd.DataFrame()
    temp3['Name'] = df.Name.tolist()
    for c in ['gene_id_known', 'gene_id_known_rt']:
        for feat in ['End', 'Start']:
            temp2 = pd.DataFrame() 
            for g in df[c].unique().tolist():
                if feat == 'Start':
                    other_feat = 'End'
                elif feat == 'End':
                    other_feat = 'Start'
                temp_nov = df.loc[df[c] == g].copy(deep=True)
                temp_nov = temp_nov[['Name', 'Chromosome', feat, 'Strand']]
                temp_nov[other_feat] = temp_nov[feat]+1 

                temp_known = known_df.loc[known_df.gene_id == g].copy(deep=True)
                temp_known = temp_known[['Chromosome', feat, 'Strand']]
                temp_known[other_feat] = temp_known[feat]+1 


                temp_nov = pr.PyRanges(temp_nov)
                temp_known = pr.PyRanges(temp_known)

                # pyranges join and get closest; concat in w/ temp2s
                temp_nov = temp_nov.k_nearest(temp_known, k=1, 
                                   overlap=True,
                                   how=None, suffix='_known').df
                try:
                    temp_nov = temp_nov[['Name', 'Distance']].drop_duplicates()
                except:
                    import pdb; pdb.set_trace()
                d_col = f'{c}_{feat}_dist'
                temp_nov.rename({'Distance': d_col}, axis=1, inplace=True)        
                # if len(temp_nov.loc[temp_nov[d_col].isnull()].index) > 0:
                #     import pdb; pdb.set_trace()
                temp2 = pd.concat([temp2, temp_nov], axis=0)
            temp3 = temp3.merge(temp2, on='Name', how='left') 
            
    df = df.merge(temp3, how='left', on='Name')  
    
    # add abs. values
    cols = ['gene_id_known_End_dist',
            'gene_id_known_Start_dist',
            'gene_id_known_rt_End_dist',
            'gene_id_known_rt_Start_dist']
    for c in cols:
        new_col = f'{c}_abs'
        df[new_col] = df[c].abs()
    
    # implement DecisionTree classifier in pandas format w/ thresholds
    df['fix'] = False

    # start of closest rt transcript must be > 510 bp away
    inds = df.loc[df.gene_id_known_rt_Start_dist_abs > 510].index
    df.loc[inds, 'fix'] = True

    # end of the closest rt transcript must be > 45 kbp away
    inds = df.loc[df.gene_id_known_rt_End_dist_abs > 45000].index
    df.loc[inds, 'fix'] = True

    # % of sss annotated in best matching non-rt transcripts
    # must be >94
    # % of sss annotated in actual transcript mush be > 84%
    inds = df.loc[(df.perc_supp_annot_ss > 94)&(df.perc_supp_annot_ss > 84)].index
    df.loc[inds, 'fix'] = True

    df[['fix_needed', 'fix', 'Name']].groupby(['fix_needed', 'fix']).count().reset_index()
    
    gtf_df = pr.read_gtf(talon_gtf, duplicate_attr=True, rename_attr=True)
    gtf_df = gtf_df.df
    gtf_df['gid_stable'] = cerberus.get_stable_gid(gtf_df, 'gene_id')

    ref_gtf_df = pr.read_gtf(ref_gtf, duplicate_attr=True, rename_attr=True)
    ref_gtf_df = ref_gtf_df.df
    ref_gtf_df['gid_stable'] = cerberus.get_stable_gid(ref_gtf_df, 'gene_id')
    
    # update the genes that need fixing to the genes that were non-rt that they 
    # intersected with
    gene_cols = ['Source', 'gene_id', 'gene_name', 'gene_status', 'gene_type', 'talon_gene', 'havana_gene', 'level',
                 'antisense_gene', 'gene_antisense_to_IDs', 'intergenic_novel', 'fusion_novel']
    talon_gene_cols = ['gene_status', 'talon_gene', 'antisense_gene', 'gene_antisense_to_IDs', 'intergenic_novel', 'fusion_novel']
    for ind, entry in temp.iterrows():
        t = entry.Name
        inds = gtf_df.loc[gtf_df.transcript_id==t].index
        g = entry.gene_id_known
        dummy_gene_entry = gtf_df.loc[(gtf_df.gid_stable==g)&(gtf_df.Feature=='gene')]
        if len(dummy_gene_entry.index) == 0:
            # have to pull from reference gtf instead and add corresponding gene entry
            dummy_gene_entry = ref_gtf_df.loc[(ref_gtf_df.gid_stable==g)&(ref_gtf_df.Feature=='gene')]
            dummy_gene_entry[talon_gene_cols] = np.nan
            gtf_df = pd.concat([gtf_df, dummy_gene_entry], axis=0)
        try:
            assert len(dummy_gene_entry.index) == 1
        except:
            import pdb; pdb.set_trace()
        for c in gene_cols:
            gtf_df.loc[inds, c] = dummy_gene_entry[c].values[0]
        
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
    pr.PyRanges(gtf_temp).to_gtf(ofile_gtf)
    
    # also update the unfiltered abundance table
    assert len(df.loc[(df.perc_supp_ss==100)&(df.perc_supp_annot_ss==100)].index) == 0
    temp = df[['Name', 'gene_id_known']]
    temp.rename({'Name': 'annot_transcript_id',
                 'gene_id_known': 'gid_stable'},
                axis=1, inplace=True)
    
    df = pd.read_csv(talon_ab, sep='\t')
    print(len(df.index))
    df['gid_stable'] = cerberus.get_stable_gid(df, 'annot_gene_id')
    df[df.columns[:11]].head()
    gene_cols = ['gene_ID', 'annot_gene_id', 'annot_gene_name', 'gene_novelty', 'gid_stable']
    gene_info = df[gene_cols].drop_duplicates()
    gene_info.head()
    temp = temp.merge(gene_info, how='left', on='gid_stable')
    temp.drop('gid_stable', axis=1, inplace=True)    
    df = df.merge(temp, how='left', on='annot_transcript_id', suffixes=('', '_fix'))
    fix_inds = df.loc[df.gene_ID_fix.notnull()].index
    gene_cols = [g for g in gene_cols if g!='gid_stable']
    for g in gene_cols:
        fix_col = f'{g}_fix'
        df.loc[fix_inds, g] = df.loc[fix_inds, fix_col]
        df.drop(fix_col, axis=1, inplace=True)
    df.to_csv(ofile_ab, sep='\t', index=False)
    
    