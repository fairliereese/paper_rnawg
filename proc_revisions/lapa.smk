
# filt_t_novs = ['Known', 'NIC', 'NNC', 'ISM_rescue']
# filt_g_novs = ['Known', 'Intergenic', 'Fusion']
# filt_spikes = True

filt_t_novs = config['params']['lapa']['filt_t_novs']
filt_g_novs = config['params']['lapa']['filt_g_novs']
filt_spikes = config['params']['lapa']['filt_spikes']

# python functions

def get_ids_from_pass_list(filt_list):
    filt_df = pd.read_csv(filt_list, sep='\t')
    gids = filt_df.gid.unique().tolist()
    tids = filt_df.tid.unique().tolist()
    return gids, tids

def ab_add_rescue_ism_cat(ab):
    """
    Update LAPA abundance w/ ISM rescue category for those
    that were assigned new ends from LAPA
    """
    df = pd.read_csv(ab, sep='\t')
    df.loc[(df.annot_transcript_id.str.contains('#'))&(df.transcript_novelty=='ISM'), 'transcript_novelty'] = 'ISM_rescue'
    return df

def filter_lapa_on_nov(df,
                        t_novs=['Known', 'NIC', 'NNC', 'ISM_rescue'],
                        g_novs=['Known']):
    """
    Filter LAPA output based on gene and transcript novelty.

    Input:
        df (pandas df): Abundance table from LAPA
        t_nov (list of str): Transcript novelty categories to keep
        g_nov (list of str): Gene novelty categories to keep

    Returns:
        filt_df (pandas df): DataFrame of gene id, transcript id passing filt
    """
    df = df.loc[df.transcript_novelty.isin(t_novs)]
    df = df.loc[df.gene_novelty.isin(g_novs)]
    filt_df = df[['annot_gene_id', 'annot_transcript_id']].drop_duplicates()
    filt_df = filt_df.rename({'annot_gene_id':'gid',
                              'annot_transcript_id': 'tid'}, axis=1)
    return filt_df

def filt_lapa_ab(ab, filt_list):
    """
    Filter LAPA abundance using a TALON-style pass list
    """
    df = pd.read_csv(ab, sep='\t')
    gids, tids = get_ids_from_pass_list(filt_list)
    df = df.loc[(df.annot_gene_id.isin(gids))&(df.annot_transcript_id.isin(tids))]
    return df

def filter_spikes(gtf):
    """
    Filter LAPA output based on SIRV / ERCC status

    Input:
        gtf (str): GTF path from LAPA
    Returns:
        filt_df (pandas df): DataFrame of gene id, transcript id passing filt
    """
    df = pr.read_gtf(gtf, as_df=True)
    df = df.loc[~df.Chromosome.str.contains('SIRV')]
    df = df.loc[~df.Chromosome.str.contains('ERCC')]
    df = df.loc[df.Feature == 'transcript']
    filt_df = df[['gene_id', 'transcript_id']].drop_duplicates()
    filt_df = filt_df.rename({'gene_id':'gid',
                              'transcript_id':'tid'}, axis=1)
    return filt_df

def filter_lapa(ab,
                gtf,
                t_novs,
                g_novs,
                filt_spikes,
                ofile):
    """
    Filter LAPA transcripts and output a pass list of
    passed transcripts that are in the list of valid
    transcript novelites, gene novelties, and depending
    on whether or not we're filtering out spike-ins.

    Parameters:
        ab (str): Path to abundance file
        gtf (str): Path to GTF
        t_novs (list of str): List of transcript novelties
            to retain transcripts from
        g_novs (list of str): List of gene novelties to
            retain transcripts from
        filt_spikes (bool): Whether to remove spikeins
        ofile (str): Path / name of output file
    """

    # filter based on novelty after defining
    # rescue ISMS
    df = ab_add_rescue_ism_cat(ab)
    filt_df = filter_lapa_on_nov(df,
                                 t_novs,
                                 g_novs)

    # filter out spike-ins
    if filt_spikes:
        temp = filter_spikes(gtf)
        filt_df = filt_df.merge(temp, how='inner')

    filt_df.to_csv(ofile, index=False, sep='\t')

def get_lapa_run_info(wc, df, cfg_entry, dataframe=False):
    """
    Get all files for a lapa run

    Parameters:
        dataframe (bool): False if it should return just the
            list of input files for the lapa run
            True if it should return the whole DF
    """
    temp = df.copy(deep=True)
    temp = temp.loc[(temp.species==wc.species)]
    datasets = temp.dataset.tolist()
    species = temp.species.tolist()
    files = expand(cfg_entry,
                   zip,
                   dataset=datasets,
                   species=species)
    temp['lapa_file'] = files
    if not dataframe:
        return files
    else:
        return temp

def filt_lapa_gtf(gtf, filt_list):
    """
    Filter LAPA GTF using a TALON-style pass list
    """
    gtf = pr.read_gtf(gtf).as_df()
    gids, tids = get_ids_from_pass_list(filt_list)

    # first filter on tids
    gtf = gtf.loc[(gtf.transcript_id.isin(tids))|(gtf.Feature=='gene')]

    # then filter on gids
    gtf = gtf.loc[(gtf.gene_id.isin(gids))]

    gtf = pr.PyRanges(gtf)
    return gtf

def get_lapa_settings(wc, lapa_ends, kind):
    """
    Get the command name or file output name
    given whether we're running LAPA in tss or tes mode

    Parameters:
        kind (str): {'temp_file', 'lapa_cmd'}
    """
    lapa_ends = expand(lapa_ends,
                       zip,
                       species=wc.species,
                       end_mode=wc.end_mode)[0]
    if kind == 'temp_file':
        temp = os.path.dirname(lapa_ends)+'/'
        if wc.end_mode == 'tes':
            temp += 'polyA_clusters.bed'
        elif wc.end_mode == 'tss':
            temp += 'tss_clusters.bed'
        return temp
    elif kind == 'lapa_cmd':
        if wc.end_mode == 'tes':
            return 'lapa'
        elif wc.end_mode == 'tss':
            return 'lapa_tss'

################################################################################
################################ Template rules ################################
################################################################################
rule lapa_link:
    resources:
        threads = 1,
        mem_gb = 256
    shell:
        """
        lapa_link_tss_to_tes \
            --alignment {input.annot} \
            --lapa_dir {params.tes_dir} \
            --lapa_tss_dir {params.tss_dir} \
            --output {output.links}
        """

rule lapa_call_ends:
    resources:
        threads = 4,
        mem_gb = 32
    shell:
        """
        rm -rf {params.opref}
        {params.lapa_cmd} \
            --alignment {input.config} \
            --fasta {input.fa} \
            --annotation {input.gtf} \
            --chrom_sizes {input.chrom_sizes} \
            --output_dir {params.opref}
        if [ {params.lapa_end_temp} != {output.ends} ]
        then
            cp {params.lapa_end_temp} {output.ends}
        fi
        """

rule lapa_correct_talon:
    resources:
        threads = 1,
        mem_gb = 128
    shell:
        """
        lapa_correct_talon \
                --links {input.links} \
                --read_annot {input.annot} \
                --gtf_input {input.gtf} \
                --gtf_output {output.gtf} \
                --abundance_input {input.filt_ab} \
                --abundance_output {output.filt_ab} \
                --keep_unsupported
        """

rule lapa_filt:
    resources:
        threads = 1,
        mem_gb = 32
    run:
        filter_lapa(input.filt_ab,
                    input.gtf,
                    params.t_novs,
                    params.g_novs,
                    params.filt_spikes,
                    output.filt_list)


rule lapa_filt_ab:
    resources:
        threads = 4,
        mem_gb = 32
    run:
        df = filt_lapa_ab(input.ab,
                          input.filt_list)
        df.to_csv(output.ab, sep='\t', index=False)

rule lapa_filt_gtf:
    resources:
        threads = 4,
        mem_gb = 32
    run:
        gtf = filt_lapa_gtf(input.gtf,
                            input.filt_list)
        gtf.to_gtf(output.gtf)



################################################################################
################################# LAPA #########################################
################################################################################
rule lapa_config:
    input:
        files = lambda wc:get_lapa_run_info(wc, lr_df,
                config['lr']['talon']['bam_sort'],
                dataframe=False)
    resources:
        threads = 1,
        mem_gb = 1
    params:
        df = lr_df
    output:
        config = config['lr']['lapa']['config']
    run:
        df = get_lapa_run_info(wildcards, params.df,
                    config['lr']['talon']['bam_sort'],
                    dataframe=True)
        # df = df[['sample', 'dataset', 'lapa_file']].copy(deep=True)
        # df.columns = ['sample', 'dataset', 'path']

        df = df[['dataset', 'sample', 'lapa_file']].copy(deep=True)

        # argggggh rename these files in the opposite manner because
        # hasan and I have opposite definitions of "sample" and "dataset"
        df.columns = ['sample', 'dataset', 'path']

        df.to_csv(output.config, sep=',', index=False)

use rule lapa_call_ends as lapa_call_ends_full with:
    input:
        config = config['lr']['lapa']['config'],
        fa = config['ref']['talon']['fa'],
        gtf = config['ref']['lapa']['gtf'],
        chrom_sizes = config['ref']['talon']['chrom_sizes']
    params:
        opref = config['lr']['lapa']['ends'].rsplit('/', maxsplit=1)[0]+'/',
        lapa_cmd = lambda wc:get_lapa_settings(wc,
                                config['lr']['lapa']['ends'],
                                'lapa_cmd'),
        lapa_end_temp = lambda wc:get_lapa_settings(wc,
                            config['lr']['lapa']['ends'],
                            'temp_file')
    output:
        ends = config['lr']['lapa']['ends']

use rule lapa_link as lapa_link_full with:
    input:
        annot = config['lr']['talon']['full_annot'],
        tss = expand(config['lr']['lapa']['ends'], end_mode='tss', allow_missing=True)[0],
        tes = expand(config['lr']['lapa']['ends'], end_mode='tes', allow_missing=True)[0]
    params:
        tss_dir = expand(config['lr']['lapa']['ends'], end_mode='tss', allow_missing=True)[0].rsplit('/', maxsplit=1)[0]+'/',
        tes_dir = expand(config['lr']['lapa']['ends'], end_mode='tes', allow_missing=True)[0].rsplit('/', maxsplit=1)[0]+'/'
    output:
        links = config['lr']['lapa']['links']

use rule lapa_correct_talon as lapa_correct_talon_full with:
    input:
        gtf = config['lr']['talon']['fusion_fix']['gtf'],
        filt_ab = config['lr']['talon']['filt_ab'],
        annot = config['lr']['talon']['full_annot'],
        links = config['lr']['lapa']['links']
    output:
        gtf = config['lr']['lapa']['gtf'],
        filt_ab = config['lr']['lapa']['filt_ab']

use rule lapa_filt as lapa_filt_full with:
    input:
        filt_ab = config['lr']['lapa']['filt_ab'],
        gtf = config['lr']['lapa']['gtf']
    params:
        t_novs = filt_t_novs,
        g_novs = filt_g_novs,
        filt_spikes = filt_spikes
    output:
        filt_list = config['lr']['lapa']['filt']['pass_list']

use rule lapa_filt_ab as lapa_filt_ab_full with:
    input:
        ab = config['lr']['lapa']['filt_ab'],
        filt_list = config['lr']['lapa']['filt']['pass_list']
    output:
        ab = config['lr']['lapa']['filt']['filt_ab']

use rule lapa_filt_gtf as lapa_filt_gtf_full with:
    input:
        gtf = config['lr']['lapa']['gtf'],
        filt_list = config['lr']['lapa']['filt']['pass_list']
    output:
        gtf = config['lr']['lapa']['filt']['gtf']
