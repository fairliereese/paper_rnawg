import cerberus
import pandas as pd
import shutil

end_modes = ['tss', 'tes']


wildcard_constraints:
    end_mode='|'.join([re.escape(x) for x in end_modes]),

################################################################################
#################### Get triplet features from GTF #############################
################################################################################
rule cerb_gtf_to_bed:
    resources:
        mem_gb = 64,
        threads = 1
    run:
        cerberus.gtf_to_bed(input.gtf,
                            wildcards.end_mode,
                            output.ends,
                            dist=params.dist,
                            slack=params.slack)

rule cerb_gtf_to_ics:
    resources:
        mem_gb = 64,
        threads = 1
    run:
        cerberus.gtf_to_ics(input.gtf,
                            output.ics)

################################################################################
######################### Cerberus annotation ##################################
################################################################################

rule cerb_write_ref:
  resources:
      threads = 4,
      mem_gb = 64
  run:
      cerberus.write_reference(input.tss,
                               input.tes,
                               input.ic,
                               output.h5)

rule cerb_annot:
   resources:
       mem_gb = 64,
       threads = 16
   run:
       # if we're dealing w/ mouse and gtex, we don't need this file
       # cause there is no gtex. just copy
       if params.source == 'gtex' and wildcards.species == 'mouse':
           shutil.copyfile(input.h5, output.h5)
       else:
           cerberus.annotate_transcriptome(input.gtf,
                                           input.h5,
                                           params.source,
                                           params.gene_source,
                                           output.h5)

rule cerb_gtf_ids:
    resources:
       mem_gb = 64,
       threads = 16
    run:
       cerberus.replace_gtf_ids(input.h5,
                                input.gtf,
                                params.source,
                                params.update_ends,
                                params.agg,
                                output.gtf)

rule cerb_ab_ids:
    resources:
        mem_gb = 64,
        threads = 16
    run:
        cerberus.replace_ab_ids(input.ab,
                                input.h5,
                                params.source,
                                params.agg,
                                output.ab)


#################################################################################
##################### Ref. ends / ics for Cerberus #############################
################################################################################

# old refs
use rule cerb_gtf_to_bed as cerb_get_gtf_ends_ref with:
    input:
        gtf = config['ref']['talon']['gtf']
    output:
        ends = config['ref']['cerberus']['ends']
    params:
        slack = lambda wc:config['params']['cerberus'][wc.end_mode]['slack'],
        dist = lambda wc:config['params']['cerberus'][wc.end_mode]['dist']

use rule cerb_gtf_to_ics as cerb_get_gtf_ics_ref with:
    input:
        gtf = config['ref']['talon']['gtf']
    output:
        ics = config['ref']['cerberus']['ics']

# new refs
use rule cerb_gtf_to_bed as cerb_get_gtf_ends_ref_new with:
    input:
        gtf = config['ref']['new_gtf']
    output:
        ends = config['ref']['cerberus']['new_ends']
    params:
        slack = lambda wc:config['params']['cerberus'][wc.end_mode]['slack'],
        dist = lambda wc:config['params']['cerberus'][wc.end_mode]['dist']

use rule cerb_gtf_to_ics as cerb_get_gtf_ics_ref_new with:
    input:
        gtf = config['ref']['new_gtf']
    output:
        ics = config['ref']['cerberus']['new_ics']

# data refs -- from lapa
use rule cerb_gtf_to_bed as cerb_get_gtf_ends_lr with:
    input:
        gtf = config['lr']['lapa']['filt']['gtf']
    output:
        ends = config['lr']['cerberus']['ends']
    params:
        slack = lambda wc:config['params']['cerberus'][wc.end_mode]['slack'],
        dist = lambda wc:config['params']['cerberus'][wc.end_mode]['dist']

use rule cerb_gtf_to_ics as cerb_get_gtf_ics_lr with:
    input:
        gtf = config['lr']['lapa']['filt']['gtf']
    output:
        ics = config['lr']['cerberus']['ics']

################################################################################
######################### Cerberus aggregation #################################
################################################################################

#### agg stuff
rule cerb_agg_human_tss_config:
    input:
        v40 = expand(config['ref']['cerberus']['new_ends'],
                     species='human',
                     end_mode='tss')[0],
        v29 = expand(config['ref']['cerberus']['ends'],
                     species='human',
                     end_mode='tss')[0],
        lapa = expand(config['lr']['cerberus']['ends'],
                      species='human',
                      end_mode='tss')[0],
        gtex = expand(config['gtex']['cerberus']['ends'],
                      end_mode='tss',
                      species='human')[0],
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
        cfg = expand(config['lr']['cerberus']['agg_ends_cfg'],
                     species='human',
                     end_mode='tss')[0]
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

rule cerb_agg_human_tes_config:
    input:
        v40 = expand(config['ref']['cerberus']['new_ends'],
                     species='human',
                     end_mode='tes')[0],
        v29 = expand(config['ref']['cerberus']['ends'],
                     species='human',
                     end_mode='tes')[0],
        lapa = expand(config['lr']['cerberus']['ends'],
                      species='human',
                      end_mode='tes')[0],
        gtex = expand(config['gtex']['cerberus']['ends'],
                      end_mode='tes',
                      species='human')[0],
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
        cfg = expand(config['lr']['cerberus']['agg_ends_cfg'],
                     species='human',
                     end_mode='tes')[0]
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

rule cerb_agg_ics_human_config:
    input:
        v40 = expand(config['ref']['cerberus']['new_ics'],
                     species='human')[0],
        v29 = expand(config['ref']['cerberus']['ics'],
                     species='human')[0],
        lapa = expand(config['lr']['cerberus']['ics'],
                     species='human')[0],
        gtex = expand(config['gtex']['cerberus']['ics'],
                      species='human')[0]
    params:
        sources = ['v40', 'v29', 'lapa', 'gtex'],
        refs = [True, True, False, False]
    resources:
        mem_gb = 1,
        threads = 1
    output:
        cfg = expand(config['lr']['cerberus']['agg_ics_cfg'],
                     species='human')[0]
    run:
        files = [input.v40, input.v29, input.lapa, input.gtex]
        refs = params.refs
        sources = params.sources
        df = pd.DataFrame()
        df['fname'] = files
        df['refs'] = refs
        df['source'] = sources
        df.to_csv(output.cfg, header=None, index=False, sep=',')

rule cerb_agg_mouse_tss_config:
    input:
        vM21 = expand(config['ref']['cerberus']['ends'],
                     species='mouse',
                     end_mode='tss')[0],
        vM25 = expand(config['ref']['cerberus']['new_ends'],
                     species='mouse',
                     end_mode='tss')[0],
        lapa = expand(config['lr']['cerberus']['ends'],
                    species='mouse',
                    end_mode='tss')[0],
        fantom_cage = expand(config['fantom']['bed'],
                             species='mouse')[0],
        pls = expand(config['ccre']['bed_format'],
                     species='mouse',
                     ccre_type='pls')[0],
        pels = expand(config['ccre']['bed_format'],
                      species='mouse',
                      ccre_type='pels')[0],
        dels = expand(config['ccre']['bed_format'],
                      species='mouse',
                      ccre_type='dels')[0],
        h3k4me3 = expand(config['ccre']['bed_format'],
                    species='mouse',
                    ccre_type='ca_h3k4me3')[0],
        pol2 = expand(config['pol2']['merged'],
                      species='mouse')[0]
    resources:
        mem_gb = 1,
        threads = 1
    params:
        add_ends = [True, True, True,
                  False, False, False, False, False,
                  False],
        refs = [True, True, False,
                  False, False, False, False, False,
                  False],
        sources = ['vM25', 'vM21', 'lapa',
                 'fantom_cage', 'pls', 'pels', 'dels', 'h3k4me3',
                 'pol2']
    output:
        cfg = expand(config['lr']['cerberus']['agg_ends_cfg'],
                   species='mouse',
                   end_mode='tss')[0]
    run:
        files = [input.vM25,
               input.vM21,
               input.lapa,
               input.fantom_cage,
               input.pls,
               input.pels,
               input.dels,
               input.h3k4me3,
               input.pol2]
        df = pd.DataFrame()
        df['fname'] = files
        df['add_ends'] = params.add_ends
        df['refs'] = params.refs
        df['sources'] = params.sources
        df.to_csv(output.cfg, sep=',', header=None, index=False)

rule cerb_agg_mouse_tes_config:
    input:
        vM21 = expand(config['ref']['cerberus']['ends'],
                     species='mouse',
                     end_mode='tes')[0],
        vM25 = expand(config['ref']['cerberus']['new_ends'],
                     species='mouse',
                     end_mode='tes')[0],
        lapa = expand(config['lr']['cerberus']['ends'],
                    species='mouse',
                    end_mode='tes')[0],
        pas = expand(config['pas']['ends_formatted'],
                species='mouse',
                end_mode='tes')[0],
        atlas = expand(config['polya_atlas']['bed_formatted'],
                       species='mouse')[0]
    resources:
        mem_gb = 1,
        threads = 1
    params:
        add_ends = [True, True, True,
                   False, False],
        refs = [True, True, False,
                   False, False],
        sources = ['vM25', 'vM21', 'lapa',
                   'pas', 'polya_atlas']
    output:
        cfg = expand(config['lr']['cerberus']['agg_ends_cfg'],
                     species='mouse',
                     end_mode='tes')[0]
    run:
        files = [input.vM25,
                input.vM21,
                input.lapa,
                input.pas,
                input.atlas]
        df = pd.DataFrame()
        df['fname'] = files
        df['add_ends'] = params.add_ends
        df['refs'] = params.refs
        df['sources'] = params.sources
        df.to_csv(output.cfg, sep=',', header=None, index=False)


rule cerb_agg_ics_mouse_config:
    input:
        vM25 = expand(config['ref']['cerberus']['new_ics'],
                     species='mouse')[0],
        vM21 = expand(config['ref']['cerberus']['ics'],
                     species='mouse')[0],
        lapa = expand(config['lr']['cerberus']['ics'],
                     species='mouse')[0]
    params:
        sources = ['vM25', 'vM21', 'lapa'],
        refs = [True, True, False]
    resources:
        mem_gb = 1,
        threads = 1
    output:
        cfg = expand(config['lr']['cerberus']['agg_ics_cfg'],
                     species='mouse')[0]
    run:
        files = [input.vM25, input.vM21, input.lapa]
        refs = params.refs
        sources = params.sources
        df = pd.DataFrame()
        df['fname'] = files
        df['refs'] = refs
        df['source'] = sources
        df.to_csv(output.cfg, header=None, index=False, sep=',')

rule cerb_agg_ends:
    input:
        cfg = config['lr']['cerberus']['agg_ends_cfg']
    resources:
        mem_gb = 64,
        threads = 1
    params:
        agg_slack = lambda wc:config['params']['cerberus'][wc.end_mode]['agg_slack']
    output:
        bed = config['lr']['cerberus']['agg_ends']
    shell:
        """
        cerberus agg_ends \
            --input {input.cfg} \
            --mode {wildcards.end_mode} \
            --slack {params.agg_slack} \
            -o {output.bed}
        """

rule cerb_agg_ics:
    input:
        cfg = config['lr']['cerberus']['agg_ics_cfg']
    resources:
        mem_gb = 28,
        threads = 1
    output:
        ics = config['lr']['cerberus']['agg_ics']
    shell:
        "cerberus agg_ics \
            --input {input.cfg} \
            -o {output.ics}"

################################################################################
######################### Cerberus annotation ##################################
################################################################################
rule cerberus_write_ref:
    input:
        tss = expand(config['lr']['cerberus']['agg_ends'],
                     end_mode='tss',
                     allow_missing=True)[0],
        tes = expand(config['lr']['cerberus']['agg_ends'],
                     end_mode='tes',
                     allow_missing=True)[0],
        ics = rules.cerb_agg_ics.output.ics
    resources:
        mem_gb = 56,
        threads = 1
    output:
        h5 = config['lr']['cerberus']['ca']
    run:
        cerberus.write_reference(input.tss,
                                 input.tes,
                                 input.ics,
                                 output.h5)

# rule cerberus_annotate:
#     resources:
#         mem_gb = 56,
#         threads = 1
#     shell:
#         """cerberus annotate_transcriptome \
#             --gtf {input.gtf} \
#             --h5 {input.ref} \
#             --source {params.source} \
#             --gene_source {params.gene_source}
#             -o {output.annot}
#         """

################################################################################
########################## Cerberus annotation #################################
################################################################################

use rule cerb_annot as cerberus_annotate_new_ref with:
    input:
        gtf = config['ref']['new_gtf'],
        h5 = config['lr']['cerberus']['ca']
    params:
        source = lambda wc:config['ref'][wc.species]['new_gtf_ver'],
        gene_source = None
    output:
        h5 = config['ref']['cerberus']['new_ca']

use rule cerb_annot as cerberus_annotate_ref with:
    input:
        gtf = config['ref']['talon']['gtf'],
        h5 = config['ref']['cerberus']['new_ca']
    params:
        source = lambda wc:config['ref'][wc.species]['gtf_ver'],
        gene_source = lambda wc:config['ref'][wc.species]['new_gtf_ver']
    output:
        h5 = config['ref']['cerberus']['ca']

use rule cerb_annot as cerberus_annotate_gtex with:
    input:
        gtf = config['gtex']['filt_gtf'],
        h5 = config['ref']['cerberus']['ca']
    params:
        source = 'gtex',
        gene_source = lambda wc:config['ref'][wc.species]['new_gtf_ver']
    output:
        h5 = config['gtex']['cerberus']['ca']

use rule cerb_annot as cerberus_annotate_lr with:
    input:
        gtf = config['lr']['lapa']['filt']['gtf'],
        h5 = config['gtex']['cerberus']['ca']
    params:
        source = 'lapa',
        gene_source = lambda wc:config['ref'][wc.species]['new_gtf_ver']
    output:
        h5 = config['lr']['cerberus']['ca_annot']

################################################################################
####################### Cerberus annotation -- mouse ###########################
################################################################################

# use rule cerb_annot as cerberus_annotate_mouse_new_ref with:
#     input:
#         gtf = expand(config['ref']['new_gtf'],
#                      species='mouse')[0],
#         ref = expand(config['lr']['cerberus']['ca'],
#                      species='mouse')[0]
#     params:
#         source = lambda wc:config['ref'][wc.species]['new_gtf_ver'],
#         gene_source = None
#     output:
#         h5 = config['ref']['cerberus']['new_ca']
#
# use rule cerb_annot as cerberus_annotate_mouse_ref with:
#     input:
#         gtf = expand(config['ref']['talon']['gtf'],
#                      species='mouse')[0],
#         ref = expand(config['ref']['cerberus']['new_ca'],
#                      species='mouse')[0]
#     params:
#         source = lambda wc:config['ref'][wc.species]['gtf_ver'],
#         gene_source = None
#     output:
#         h5 = config['ref']['cerberus']['ca']
#
# use rule cerb_annot as cerberus_annotate_mouse_lr with:
#     input:
#         gtf = expand(config['lr']['lapa']['filt']['gtf'],
#                      species='mouse')[0],
#         ref = expand(config['ref']['cerberus']['ca'],
#                      species='mouse')[0]
#     params:
#         source = 'lapa',
#         gene_source = lambda wc:config['ref'][wc.species]['gtf_ver']
#     output:
#         h5 = config['lr']['cerberus']['ca_annot']

################################################################################
###################### Cerberus update abundance ###############################
################################################################################

use rule cerb_ab_ids as cerb_ab_ids_lr with:
    input:
        h5 = config['lr']['cerberus']['ca_annot'],
        ab = config['lr']['lapa']['filt']['filt_ab']
    params:
        source = 'lapa',
        agg = True
    output:
        ab = config['lr']['cerberus']['ab']

use rule cerb_gtf_ids as cerb_gtf_ids_lr with:
    input:
        h5 = config['lr']['cerberus']['ca_annot'],
        gtf = config['lr']['lapa']['filt']['gtf']
    params:
        source = 'lapa',
        update_ends = True,
        agg = True
    output:
        gtf = config['lr']['cerberus']['gtf']

use rule cerb_gtf_ids as cerb_gtf_ids_ref with:
    input:
        h5 = config['lr']['cerberus']['ca_annot'],
        gtf = config['ref']['talon']['gtf']
    params:
        source = lambda wc:config['ref'][wc.species]['gtf_ver'],
        update_ends = True,
        agg = True
    output:
        gtf = config['ref']['cerberus']['gtf']

use rule cerb_gtf_ids as cerb_gtf_ids_new_ref with:
    input:
        h5 = config['lr']['cerberus']['ca_annot'],
        gtf = config['ref']['new_gtf']
    params:
        source = lambda wc:config['ref'][wc.species]['new_gtf_ver'],
        update_ends = True,
        agg = True
    output:
        gtf = config['ref']['cerberus']['new_gtf']


################################################################################
######################### Get ICs from TALON GTF ###############################
################################################################################

use rule cerb_gtf_to_ics as cerb_get_gtf_ics_talon with:
    input:
        gtf = config['lr']['talon']['fusion_fix']['gtf']
    output:
        ics = config['lr']['talon']['ics']


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
        ab = config['lr']['cerberus']['ab'],
        ca = config['lr']['cerberus']['ca_annot']
    resources:
        threads = 1,
        mem_gb = 32
    output:
        filt_ab = config['lr']['cerberus']['filt_ab']
    run:
        filt_unsup_ism(input.ab, input.ca, wildcards, output.filt_ab)

################################################################################
# get ends for milad's tss prediction

def get_cerb_tss(cerberus_h5, filt_ab, species, params):
    datasets = get_datasets(species=species)
    ab_df = pd.read_csv(filt_ab, sep='\t')
    df = get_det_table(ab_df,
                       how='tss',
                       min_tpm=params.min_tpm,
                       groupby='library',
                       gene_subset=None)
    df = df.transpose()
    tpm_df, _ = get_tpm_table(ab_df,
                   how='tss',
                   min_tpm=params.min_tpm,
                   groupby='library',
                   gene_subset=None)
    ca = cerberus.read(cerberus_h5)
    tss_df = ca.tss.copy(deep=True)
    for d in datasets:
        temp = df.loc[df[d]==True].copy(deep=True)[[d]]
        beep = tpm_df[[d]]
        beep.rename({d: 'tpm'}, axis=1, inplace=True)
        temp = temp.merge(beep, left_index=True, right_index=True, how='left')
        temp = temp.merge(tss_df, how='left', left_index=True, right_on='Name')
        temp['dataset'] = d
        cols = ['Chromosome', 'Start', 'End', 'Name',
               'Strand', 'gene_id', 'tpm',
               'source', 'novelty', 'dataset']
        temp = temp[cols]
        # fname = f'{d}_cerberus.bed'
        fname = f'{params.opref}/{d}_cerberus_tss.bed'
        temp.to_csv(fname, sep='\t', index=False)
    return

def get_output_cerb_get_human_tss_ends(species, df):
    temp = df.loc[df.species==species]
    files = expand(config['lr']['cerberus']['dataset']['tss'],
                 zip,
                 species=temp['species'].tolist(),
                 dataset=temp['dataset'].tolist())
    return files

rule cerb_get_human_tss_ends:
    input:
        h5 = expand(config['lr']['cerberus']['ca_annot'],
                    species='human')[0],
        filt_ab = expand(config['lr']['cerberus']['filt_ab'],
                         species='human')[0]
    resources:
        threads = 1,
        mem_gb = 32
    params:
        min_tpm = 1,
        opref = expand(config['lr']['cerberus']['dataset']['tss'],
                       species='human', allow_missing=True)[0].rsplit('/', maxsplit=1)[0]
    output:
        bed = get_output_cerb_get_human_tss_ends('human', lr_df)
    run:
        get_cerb_tss(input.h5,
                     input.filt_ab,
                     'human',
                     params)

rule all_cerberus:
    input:
        expand(config['lr']['talon']['ics'], species=species),
        get_output_cerb_get_human_tss_ends('human', lr_df),
        expand(rules.cerb_gtf_ids_lr.output, species=species)


#         expand(rules.cerb_gtf_ids_new_ref.output, species=species),
#         expand(rules.cerb_gtf_ids_ref.output, species=species),

        # expand(rules.cerb_gtf_ids_new_ref.output,
        #        species=species),
        # expand(rules.cerb_gtf_ids_ref.output,
        #        species=species),
        # expand(rules.cerb_gtf_ids_lr.output,
        #       species=species)


        # expand(config['lr']['cerberus']['agg_ends'],
        #        species=species,
        #        end_mode=end_modes),
        # expand(config['lr']['cerberus']['agg_ics'],
        #       species=species)
        # expand(config['ref']['cerberus']['ends'],
        #       species=species,
        #       end_mode=end_modes),
        # expand(config['ref']['cerberus']['new_ends'],
        #     species=species,
        #     end_mode=end_modes),
        # expand(config['ref']['cerberus']['ics'],
        #       species=species),
        # expand(config['ref']['cerberus']['new_ics'],
        #     species=species),
        # expand(config['ref']['cerberus']['ends'],
        #       species=species,
        #       end_mode=end_modes),
        # expand(config['lr']['cerberus']['ics'],
        #     species=species)
