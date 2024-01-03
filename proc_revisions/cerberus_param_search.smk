import cerberus
import pandas as pd
import shutil

end_modes = ['tss', 'tes']


wildcard_constraints:
    end_mode='|'.join([re.escape(x) for x in end_modes]),

#################################################################################
##################### Ref. ends / ics for Cerberus #############################
################################################################################
def get_slack(wc):
    if wc.end_mode == 'tss':
        return wc.tss_slack
    elif wc.end_mode == 'tes':
        return wc.tes_slack

def get_dist(wc):
    if wc.end_mode == 'tss':
        return wc.tss_dist
    elif wc.end_mode == 'tes':
        return wc.tes_dist

def get_agg_dist(wc):
    if wc.end_mode == 'tss':
        return wc.tss_agg_dist
    elif wc.end_mode == 'tes':
        return wc.tes_agg_dist

# old refs
use rule cerb_gtf_to_bed as param_cerb_get_gtf_ends_ref_old with:
    input:
        gtf = config['ref']['talon']['gtf']
    output:
        ends = config['ref']['param_search']['cerberus']['ends']
    params:
        slack = lambda wc:get_slack(wc),
        dist = lambda wc:get_dist(wc)

use rule cerb_gtf_to_ics as param_cerb_get_gtf_ics_ref_old with:
    input:
        gtf = config['ref']['talon']['gtf']
    output:
        ics = config['ref']['param_search']['cerberus']['ics']

# new refs
use rule cerb_gtf_to_bed as param_cerb_get_gtf_ends_ref_new with:
    input:
        gtf = config['ref']['new_gtf']
    output:
        ends = config['ref']['param_search']['cerberus']['new_ends']
    params:
        slack = lambda wc:get_slack(wc),
        dist = lambda wc:get_dist(wc)

use rule cerb_gtf_to_ics as param_cerb_get_gtf_ics_ref_new with:
    input:
        gtf = config['ref']['new_gtf']
    output:
        ics = config['ref']['param_search']['cerberus']['new_ics']

# data refs -- from lapa
use rule cerb_gtf_to_bed as param_cerb_get_gtf_ends_lr_ref with:
    input:
        gtf = config['lr']['lapa']['filt']['gtf']
    output:
        ends = config['lr']['param_search']['cerberus']['ends']
    params:
        slack = lambda wc:get_slack(wc),
        dist = lambda wc:get_dist(wc)

use rule cerb_gtf_to_ics as param_cerb_get_gtf_ics_lr with:
    input:
        gtf = config['lr']['lapa']['filt']['gtf']
    output:
        ics = config['lr']['param_search']['cerberus']['ics']

################################################################################
######################### Cerberus aggregation #################################
################################################################################

#### agg stuff
rule param_cerb_agg_human_tss_config:
    input:
        v40 = lambda wc: expand(config['ref']['param_search']['cerberus']['new_ends'],
                     species='human',
                     end_mode='tss',
                     tss_dist=wc.tss_dist,
                     tes_dist=wc.tes_dist,
                     tss_slack=wc.tss_slack,
                     tes_slack=wc.tes_slack,
                     tss_agg_dist=wc.tss_agg_dist,
                     tes_agg_dist=wc.tes_agg_dist)[0],
        v29 = lambda wc: expand(config['ref']['param_search']['cerberus']['ends'],
                     species='human',
                     end_mode='tss',
                     tss_dist=wc.tss_dist,
                     tes_dist=wc.tes_dist,
                     tss_slack=wc.tss_slack,
                     tes_slack=wc.tes_slack,
                     tss_agg_dist=wc.tss_agg_dist,
                     tes_agg_dist=wc.tes_agg_dist)[0],
        lapa = lambda wc: expand(config['lr']['param_search']['cerberus']['ends'],
                      species='human',
                      end_mode='tss',
                      tss_dist=wc.tss_dist,
                      tes_dist=wc.tes_dist,
                      tss_slack=wc.tss_slack,
                      tes_slack=wc.tes_slack,
                      tss_agg_dist=wc.tss_agg_dist,
                      tes_agg_dist=wc.tes_agg_dist)[0],
        gtex = lambda wc: expand(config['gtex']['param_search']['cerberus']['ends'],
                      end_mode='tss',
                      species='human',
                      tss_dist=wc.tss_dist,
                      tes_dist=wc.tes_dist,
                      tss_slack=wc.tss_slack,
                      tes_slack=wc.tes_slack,
                      tss_agg_dist=wc.tss_agg_dist,
                      tes_agg_dist=wc.tes_agg_dist)[0],
        encode_cage = expand(config['cage']['merged'],
                             species='human'),
        fantom_cage = expand(config['fantom']['bed'],
                             species='human')[0],
        encode_rampage = expand(config['rampage']['merged'],
                                species='human')[0],
        pls = expand(config['ccre']['bed_format'],
                     species='human',
                     ccre_type='pls')[0],
        pels = expand(config['ccre']['bed_format'],
                      species='human',
                      ccre_type='pels')[0],
        dels = expand(config['ccre']['bed_format'],
                      species='human',
                      ccre_type='dels')[0],
        h3k4me3 = expand(config['ccre']['bed_format'],
                         species='human',
                         ccre_type='ca_h3k4me3')[0],
        lrgasp_cage = expand(config['lrgasp_cage']['merged'],
                             species='human')[0],
        encode_procap = expand(config['procap']['merged'],
                               species='human')[0],
        pol2 = expand(config['pol2']['merged'],
                      species='human')[0]
    params:
        add_ends = [True, True, True, True,
                    False, False, False,
                    False, False, False, False,
                    False, False, False],
        refs = [True, True, False, False,
                False, False, False,
                False, False, False, False,
                False, False, False],
        sources = ['v40', 'v29', 'lapa', 'gtex',
                   'encode_cage', 'fantom_cage',
                   'encode_rampage',
                   'pls', 'pels', 'dels', 'ca_h3k4me3',
                   'lrgasp_cage', 'encode_procap', 'pol2']
    resources:
        threads = 1,
        mem_gb = 1
    output:
        cfg = config['lr']['param_search']['cerberus']['agg_ends_cfg']
    run:
        files = [input.v40,
                 input.v29,
                 input.lapa,
                 input.gtex,
                 input.encode_cage,
                 input.fantom_cage,
                 input.encode_rampage,
                 input.pls,
                 input.pels,
                 input.dels,
                 input.h3k4me3,
                 input.lrgasp_cage,
                 input.encode_procap,
                 input.pol2]
        df = pd.DataFrame()
        df['fname'] = files
        df['add_ends'] = params.add_ends
        df['refs'] = params.refs
        df['sources'] = params.sources
        df.to_csv(output.cfg, sep=',', header=None, index=False)

rule param_cerb_agg_human_tes_config:
    input:
        v40 = lambda wc: expand(config['ref']['param_search']['cerberus']['new_ends'],
                     species='human',
                     end_mode='tes',
                     tss_dist=wc.tss_dist,
                     tes_dist=wc.tes_dist,
                     tss_slack=wc.tss_slack,
                     tes_slack=wc.tes_slack,
                     tss_agg_dist=wc.tss_agg_dist,
                     tes_agg_dist=wc.tes_agg_dist)[0],
        v29 = lambda wc: expand(config['ref']['param_search']['cerberus']['ends'],
                     species='human',
                     end_mode='tes',
                     tss_dist=wc.tss_dist,
                     tes_dist=wc.tes_dist,
                     tss_slack=wc.tss_slack,
                     tes_slack=wc.tes_slack,
                     tss_agg_dist=wc.tss_agg_dist,
                     tes_agg_dist=wc.tes_agg_dist)[0],
        lapa = lambda wc: expand(config['lr']['param_search']['cerberus']['ends'],
                      species='human',
                      end_mode='tes',
                      tss_dist=wc.tss_dist,
                      tes_dist=wc.tes_dist,
                      tss_slack=wc.tss_slack,
                      tes_slack=wc.tes_slack,
                      tss_agg_dist=wc.tss_agg_dist,
                      tes_agg_dist=wc.tes_agg_dist)[0],
        gtex = lambda wc: expand(config['gtex']['param_search']['cerberus']['ends'],
                      end_mode='tes',
                      species='human',
                      tss_dist=wc.tss_dist,
                      tes_dist=wc.tes_dist,
                      tss_slack=wc.tss_slack,
                      tes_slack=wc.tes_slack,
                      tss_agg_dist=wc.tss_agg_dist,
                      tes_agg_dist=wc.tes_agg_dist)[0],
        pas = expand(config['pas']['ends_formatted'],
                species='human',
                end_mode='tes')[0],
        atlas = expand(config['polya_atlas']['bed_formatted'],
                       species='human')[0]
    resources:
        mem_gb = 1,
        threads = 1
    params:
        add_ends = [True, True, True, True,
                    False, False],
        refs = [True, True, False, False,
                False, False],
        sources = ['v40', 'v29', 'lapa', 'gtex',
                   'pas', 'polya_atlas']
    output:
        cfg = config['lr']['param_search']['cerberus']['agg_ends_cfg'],
    run:
         files = [input.v40,
                  input.v29,
                  input.lapa,
                  input.gtex,
                  input.pas,
                  input.atlas]
         df = pd.DataFrame()
         df['fname'] = files
         df['add_ends'] = params.add_ends
         df['refs'] = params.refs
         df['sources'] = params.sources
         df.to_csv(output.cfg, sep=',', header=None, index=False)

rule param_cerb_agg_ics_human_config:
    input:
        v40 = lambda wc: expand(config['ref']['param_search']['cerberus']['new_ics'],
                     species='human',
                     tss_dist=wc.tss_dist,
                     tes_dist=wc.tes_dist,
                     tss_slack=wc.tss_slack,
                     tes_slack=wc.tes_slack,
                     tss_agg_dist=wc.tss_agg_dist,
                     tes_agg_dist=wc.tes_agg_dist)[0],
        v29 = lambda wc: expand(config['ref']['param_search']['cerberus']['ics'],
                     species='human',
                     tss_dist=wc.tss_dist,
                     tes_dist=wc.tes_dist,
                     tss_slack=wc.tss_slack,
                     tes_slack=wc.tes_slack,
                     tss_agg_dist=wc.tss_agg_dist,
                     tes_agg_dist=wc.tes_agg_dist)[0],
        lapa = lambda wc: expand(config['lr']['param_search']['cerberus']['ics'],
                     species='human',
                     tss_dist=wc.tss_dist,
                     tes_dist=wc.tes_dist,
                     tss_slack=wc.tss_slack,
                     tes_slack=wc.tes_slack,
                     tss_agg_dist=wc.tss_agg_dist,
                     tes_agg_dist=wc.tes_agg_dist)[0],
        gtex = lambda wc: expand(config['gtex']['param_search']['cerberus']['ics'],
                      species='human')[0]
    params:
        sources = ['v40', 'v29', 'lapa', 'gtex'],
        refs = [True, True, False, False]
    resources:
        mem_gb = 1,
        threads = 1
    output:
        cfg = config['lr']['param_search']['cerberus']['agg_ics_cfg']
    run:
        files = [input.v40, input.v29, input.lapa, input.gtex]
        refs = params.refs
        sources = params.sources
        df = pd.DataFrame()
        df['fname'] = files
        df['refs'] = refs
        df['source'] = sources
        df.to_csv(output.cfg, header=None, index=False, sep=',')

rule param_cerb_agg_ends:
    input:
        cfg = config['lr']['param_search']['cerberus']['agg_ends_cfg']
    resources:
        mem_gb = 64,
        threads = 1
    params:
        agg_slack = lambda wc:get_agg_dist(wc)
    output:
        bed = config['lr']['param_search']['cerberus']['agg_ends']
    shell:
        """
        cerberus agg_ends \
            --input {input.cfg} \
            --mode {wildcards.end_mode} \
            --slack {params.agg_slack} \
            -o {output.bed}
        """

rule param_cerb_agg_ics:
    input:
        cfg = config['lr']['param_search']['cerberus']['agg_ics_cfg']
    resources:
        mem_gb = 28,
        threads = 1
    output:
        ics = config['lr']['param_search']['cerberus']['agg_ics']
    shell:
        "cerberus agg_ics \
            --input {input.cfg} \
            -o {output.ics}"

################################################################################
########################## Cerberus annotation #################################
################################################################################

use rule cerb_annot as param_cerberus_annotate with:
    input:
        gtf = config['ref']['new_gtf'],
        h5 = config['lr']['param_search']['cerberus']['ca']
    params:
        source = lambda wc:config['ref'][wc.species]['new_gtf_ver'],
        gene_source = None
    output:
        h5 = config['ref']['param_search']['cerberus']['new_ca']

use rule cerb_annot as param_cerberus_annotate with:
    input:
        gtf = config['ref']['talon']['gtf'],
        h5 = config['ref']['param_search']['cerberus']['new_ca']
    params:
        source = lambda wc:config['ref'][wc.species]['gtf_ver'],
        gene_source = lambda wc:config['ref'][wc.species]['new_gtf_ver']
    output:
        h5 = config['ref']['param_search']['cerberus']['ca']

use rule cerb_annot as param_cerberus_annotate_gtex with:
    input:
        gtf = config['gtex']['filt_gtf'],
        h5 = config['ref']['param_search']['cerberus']['ca']
    params:
        source = 'gtex',
        gene_source = lambda wc:config['ref'][wc.species]['new_gtf_ver']
    output:
        h5 = config['gtex']['param_search']['cerberus']['ca']

use rule cerb_annot as param_cerberus_annotate_lr with:
    input:
        gtf = config['lr']['lapa']['filt']['gtf'],
        h5 = config['gtex']['param_search']['cerberus']['ca']
    params:
        source = 'lapa',
        gene_source = lambda wc:config['ref'][wc.species]['new_gtf_ver']
    output:
        h5 = config['lr']['param_search']['cerberus']['ca_annot']

################################################################################
###################### Cerberus update abundance ###############################
################################################################################

use rule cerb_ab_ids as param_cerb_ab_ids_lr with:
    input:
        h5 = config['lr']['param_search']['cerberus']['ca_annot'],
        ab = config['lr']['lapa']['filt']['filt_ab']
    params:
        source = 'lapa',
        agg = True
    output:
        ab = config['lr']['param_search']['cerberus']['ab']

use rule cerb_gtf_ids as param_cerb_gtf_ids_lr with:
    input:
        h5 = config['lr']['param_search']['cerberus']['ca_annot'],
        gtf = config['lr']['lapa']['filt']['gtf']
    params:
        source = 'lapa',
        update_ends = True,
        agg = True
    output:
        gtf = config['lr']['param_search']['cerberus']['gtf']

use rule cerb_gtf_ids as param_cerb_gtf_ids with:
    input:
        h5 = config['lr']['param_search']['cerberus']['ca_annot'],
        gtf = config['ref']['talon']['gtf']
    params:
        source = lambda wc:config['ref'][wc.species]['gtf_ver'],
        update_ends = True,
        agg = True
    output:
        gtf = config['ref']['param_search']['cerberus']['gtf']

use rule cerb_gtf_ids as param_cerb_gtf_ids with:
    input:
        h5 = config['lr']['param_search']['cerberus']['ca_annot'],
        gtf = config['ref']['new_gtf']
    params:
        source = lambda wc:config['ref'][wc.species]['new_gtf_ver'],
        update_ends = True,
        agg = True
    output:
        gtf = config['ref']['param_search']['cerberus']['new_gtf']

#
# ################################################################################
# ######################### Get ICs from TALON GTF ###############################
# ################################################################################
#
# use rule cerb_gtf_to_ics as param_cerb_get_gtf_ics_talon with:
#     input:
#         gtf = config['lr']['talon']['fusion_fix']['gtf']
#     output:
#         ics = config['lr']['talon']['ics']


################################################################################
############################# More filtering ###################################
################################################################################
def filt_unsup_ism(filt_ab, cerberus_h5, wildcards, ofile):
    species=wildcards['species']
    feat = 'tss'
    if species == 'human':
        ref_sources = ['v29', 'v40']
        support_sources = ['encode_cage', 'fantom_cage',
                           'encode_rampage', 'gtex', 'pls',
                           'encode_procap', 'lrgasp_cage', 'pol2', 'ca_h3k4me3']
    elif species == 'mouse':
        ref_sources = ['vM21', 'vM25']
        support_sources = ['h3k4me3', 'fantom_cage', 'pls', 'pol2']

    tss_df = get_feat_support(filt_ab,
                              cerberus_h5,
                              feat,
                              ref_sources,
                              support_sources,
                              min_tpm=0,
                              how=feat,
                              species=species)
    feat = 'tes'
    if species == 'human':
        support_sources = ['gtex', 'pas', 'polya_atlas']
    elif species == 'mouse':
        support_sources = ['pas', 'polya_atlas']

    tes_df = get_feat_support(filt_ab,
                            cerberus_h5,
                            feat,
                            ref_sources,
                            support_sources,
                            min_tpm=0,
                            how=feat,
                            species=species)

    df = pd.read_csv(filt_ab, sep='\t')
    df = add_feat(df, 'annot_transcript_id', 'tss')
    df = add_feat(df, 'annot_transcript_id', 'tes')
    df = add_feat(df, 'annot_transcript_id', 'ic')
    ca = cerberus.read(cerberus_h5)
    temp_ic = ca.ic.drop('ic', axis=1).rename({'Name': 'ic'}, axis=1)
    df = df.merge(temp_ic, on='ic', how='left')
    rm_tids = []
    rm_tids += df.loc[df.novelty=='Unspliced'].annot_transcript_id.tolist()
    tss_df = tss_df.rename({'Name': 'tss', 'support':'tss_support'}, axis=1)
    tes_df = tes_df.rename({'Name': 'tes', 'support':'tes_support'}, axis=1)
    df = df.merge(tss_df, how='left', on='tss')
    df = df.merge(tes_df, how='left', on='tes')

    # unsupported at both
    rm_tids += df.loc[(df.novelty=='ISM')&\
                      (df.tss_support=='Novel')&\
                      (df.tes_support=='Novel')].annot_transcript_id.tolist()
    # unsupported at tss
    rm_tids += df.loc[(df.novelty=='ISM')&\
                    (df.tss_support=='Novel')].annot_transcript_id.tolist()
    # unsupported at tes
    rm_tids += df.loc[(df.novelty=='ISM')&\
                      (df.tes_support=='Novel')].annot_transcript_id.tolist()
    keep_tids = df.loc[~df.annot_transcript_id.isin(rm_tids)].annot_transcript_id.tolist()

    # filter the abundance file
    df = pd.read_csv(filt_ab, sep='\t')
    df = df.loc[df.annot_transcript_id.isin(keep_tids)]
    df.to_csv(ofile, sep='\t', index=False)

rule cerb_filt_unsup_ism:
    input:
        ab = config['lr']['param_search']['cerberus']['ab'],
        ca = config['lr']['param_search']['cerberus']['ca_annot']
    resources:
        threads = 1,
        mem_gb = 32
    output:
        filt_ab = config['lr']['param_search']['cerberus']['filt_ab']
    run:
        filt_unsup_ism(input.ab, input.ca, wildcards, output.filt_ab)


################################################################################
################################# Swan #########################################
################################################################################
import swan_vis as swan

def make_sg(input, params, wildcards):

    # initialize
    sg = swan.SwanGraph()
    sg.add_annotation(input.annot)
    sg.add_transcriptome(input.gtf, include_isms=True)
    sg.save_graph(params.prefix)

    sg.add_abundance(input.ab)
    sg.add_abundance(input.gene_ab, how='gene')
    sg.save_graph(params.prefix)

    # add metadata and add colors
    sg.add_metadata(input.meta)
    c_dict, order = get_biosample_colors(wildcards.species)
    sg.set_metadata_colors('sample', c_dict)

    # human only settings
    if wildcards.species == 'human':
        c_dict, order = get_ad_colors()
        sg.set_metadata_colors('health_status', c_dict)
    # save
    sg.save_graph(params.prefix)

# rule swan_gene_ab_add_stable_gid:
#     input:
#         ab = config['lr']['talon']['fusion_fix']['ab']
#     resources:
#         mem_gb = 24,
#         threads = 1
#     output:
#         ab = config['lr']['talon']['ab_stable_gid']
#     run:
#         df = pd.read_csv(input.ab, sep='\t')
#         df['gid_stable'] = cerberus.get_stable_gid(df, 'annot_gene_id')
#         df['annot_gene_id'] = df['gid_stable']
#         df.drop('gid_stable', axis=1, inplace=True)
#         df.to_csv(output.ab, sep='\t', index=False)

rule swan_init:
    input:
        annot = config['ref']['param_search']['cerberus']['new_gtf'],
        ab = config['lr']['param_search']['cerberus']['filt_ab'],
        gene_ab = config['lr']['talon']['ab_stable_gid'],
        gtf = config['lr']['param_search']['cerberus']['gtf'],
        meta = config['lr']['meta']
    params:
        prefix = config['lr']['param_search']['swan']['sg'].replace('.p', '')
    resources:
        mem_gb = 64,
        threads = 1
    output:
        sg = config['lr']['param_search']['swan']['sg']
    run:
        make_sg(input, params, wildcards)
#
# ################################################################################
# # get ends for milad's tss prediction
#
# def get_cerb_tss(cerberus_h5, filt_ab, species, params):
#     datasets = get_datasets(species=species)
#     ab_df = pd.read_csv(filt_ab, sep='\t')
#     df = get_det_table(ab_df,
#                        how='tss',
#                        min_tpm=params.min_tpm,
#                        groupby='library',
#                        gene_subset=None)
#     df = df.transpose()
#     tpm_df, _ = get_tpm_table(ab_df,
#                    how='tss',
#                    min_tpm=params.min_tpm,
#                    groupby='library',
#                    gene_subset=None)
#     ca = cerberus.read(cerberus_h5)
#     tss_df = ca.tss.copy(deep=True)
#     for d in datasets:
#         temp = df.loc[df[d]==True].copy(deep=True)[[d]]
#         beep = tpm_df[[d]]
#         beep.rename({d: 'tpm'}, axis=1, inplace=True)
#         temp = temp.merge(beep, left_index=True, right_index=True, how='left')
#         temp = temp.merge(tss_df, how='left', left_index=True, right_on='Name')
#         temp['dataset'] = d
#         cols = ['Chromosome', 'Start', 'End', 'Name',
#                'Strand', 'gene_id', 'tpm',
#                'source', 'novelty', 'dataset']
#         temp = temp[cols]
#         # fname = f'{d}_cerberus.bed'
#         fname = f'{params.opref}/{d}_cerberus_tss.bed'
#         temp.to_csv(fname, sep='\t', index=False)
#     return
#
# def get_output_cerb_get_human_tss_ends(species, df):
#     temp = df.loc[df.species==species]
#     files = expand(config['lr']['param_search']['cerberus']['dataset']['tss'],
#                  zip,
#                  species=temp['species'].tolist(),
#                  dataset=temp['dataset'].tolist())
#     return files
#
# rule cerb_get_human_tss_ends:
#     input:
#         h5 = expand(config['lr']['param_search']['cerberus']['ca_annot'],
#                     species='human')[0],
#         filt_ab = expand(config['lr']['param_search']['cerberus']['filt_ab'],
#                          species='human')[0]
#     resources:
#         threads = 1,
#         mem_gb = 32
#     params:
#         min_tpm = 1,
#         opref = expand(config['lr']['param_search']['cerberus']['dataset']['tss'],
#                        species='human', allow_missing=True)[0].rsplit('/', maxsplit=1)[0]
#     output:
#         bed = get_output_cerb_get_human_tss_ends('human', lr_df)
#     run:
#         get_cerb_tss(input.h5,
#                      input.filt_ab,
#                      'human',
#                      params)

################################################################################
#################################### From analysis.smk #########################
################################################################################
rule major_isos:
    input:
        sg = rules.swan_init.output.sg,
        filt_ab = rules.swan_init.input.ab,
        g_info = config['ref']['cerberus']['new_gtf_g_info']
    resources:
        mem_gb = 16,
        threads = 8
    params:
        min_tpm = 1,
        gene_subset = 'polya'
    output:
        ofile = config['lr']['param_search']['analysis']['major_isos']
    run:
        get_major_isos(input.sg,
                       input.filt_ab,
                       wildcards.obs_col,
                       wildcards.species,
                       output.ofile,
                       min_tpm=params.min_tpm,
                       gene_subset=params.gene_subset)


rule calc_triplets:
    input:
        swan_file = config['lr']['swan']['sg'],
        h5 = config['lr']['param_search']['cerberus']['ca_annot'],
        filt_ab = config['lr']['param_search']['cerberus']['filt_ab'],
        major_isos = config['lr']['param_search']['analysis']['major_isos'],
    params:
        min_tpm = 1,
        gene_subset = 'polya',
        obs_col = 'sample'
    resources:
        threads = 1,
        mem_gb = 64
    output:
        trips = config['lr']['param_search']['cerberus']['ca_triplets'],
        tsv = config['lr']['param_search']['analysis']['triplets']
    run:
        if wildcards.species == 'human':
            calculate_human_triplets(input.swan_file,
                                input.h5,
                                input.filt_ab,
                                input.major_isos,
                                output.trips,
                                output.tsv,
                                obs_col=params.obs_col,
                                min_tpm=params.min_tpm,
                                gene_subset=params.gene_subset)
        elif wildcards.species == 'mouse':
            calculate_mouse_triplets(input.swan_file,
                              input.h5,
                              input.filt_ab,
                              input.major_isos,
                              output.trips,
                              output.tsv,
                              obs_col=params.obs_col,
                              min_tpm=params.min_tpm,
                              gene_subset=params.gene_subset)

# def get_fusion_sample_t_coords(ab, gtf, min_tpm, sample, species, ofile):
#     """
#     Get genomic start / stop coords for transcripts expressed in
#     a given sample
#     """
#
#     df = pd.read_csv(ab, sep='\t')
#     tids = df.loc[df.gene_novelty=='Fusion', 'annot_transcript_id'].tolist()
#     print(len(tids))
#     df = get_det_table(df,
#                      groupby='sample',
#                      how='iso',
#                      min_tpm=min_tpm,
#                      species=species)
#     df = df.transpose()
#     tids2 = df.loc[df[sample]==True].index.tolist()
#     tids = list(set(tids)&set(tids2))
#
#     gtf_df = pr.read_gtf(gtf).df
#     gtf_df = gtf_df.loc[gtf_df.transcript_id.isin(tids)]
#
#     gtf_df = gtf_df.loc[gtf_df.Feature=='transcript']
#     gtf_df = gtf_df[['Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'transcript_id',]]
#
#     gtf_df.to_csv(ofile, sep='\t', index=False)

# rule get_fusion_sample_t_coords:
#     input:
#         filt_ab = config['lr']['param_search']['cerberus']['filt_ab'],
#         gtf = config['lr']['param_search']['cerberus']['gtf'],
#         g_info = config['ref']['cerberus']['new_gtf_g_info']
#     resources:
#         mem_gb = 16,
#         threads = 1
#     params:
#         min_tpm = 1
#     output:
#         tsv = config['lr']['param_search']['analysis']['fusion_t_coords']
#     run:
#         get_fusion_sample_t_coords(input.filt_ab,
#             input.gtf,
#             params.min_tpm,
#             wildcards.sample,
#             wildcards.species,
#             output.tsv)

tss_dists = config['params']['param_search']['cerberus']['tss']['dist']
tes_dists = config['params']['param_search']['cerberus']['tes']['dist']
tss_slacks = config['params']['param_search']['cerberus']['tss']['slack']
tes_slacks = config['params']['param_search']['cerberus']['tes']['slack']
tss_agg_dists = config['params']['param_search']['cerberus']['tss']['agg_slack']
tes_agg_dists = config['params']['param_search']['cerberus']['tes']['agg_slack']

rule all_cerberus_param_search:
    input:
        expand(config['lr']['param_search']['analysis']['major_isos'],
               species='human',
               obs_col='sample',
               tss_dist=tss_dists,
               tes_dist=tes_dists,
               tss_slack=tss_slacks,
               tes_slack=tes_slacks,
               tss_agg_dist=tss_agg_dists,
               tes_agg_dist=tes_agg_dists),
        expand(config['lr']['param_search']['cerberus']['ca_triplets'],
               species='human',
               obs_col='sample',
               tss_dist=tss_dists,
               tes_dist=tes_dists,
               tss_slack=tss_slacks,
               tes_slack=tes_slacks,
               tss_agg_dist=tss_agg_dists,
               tes_agg_dist=tes_agg_dists)

# rule all_cerberus:
#     input:
        # expand(config['lr']['talon']['ics'], species=species),
        # get_output_cerb_get_human_tss_ends('human', lr_df)
