import cerberus
import pandas as pd

end_modes = ['tss', 'tes']
species = ['human', 'mouse']

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
        gtf = config['ref']['gtf']
    output:
        ends = config['ref']['cerberus']['ends']
    params:
        slack = lambda wc:config['params']['cerberus'][wc.end_mode]['slack'],
        dist = lambda wc:config['params']['cerberus'][wc.end_mode]['dist']

use rule cerb_gtf_to_ics as cerb_get_gtf_ics_ref with:
    input:
        gtf = config['ref']['gtf']
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
                      end_mode='tss')[0],
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
                               species='human')[0]
    params:
        add_ends = [True, True, True, True,
                    False, False, False,
                    False, False, False, False,
                    False, False],
        refs = [True, True, False, False,
                False, False, False,
                False, False, False, False,
                False, False],
        sources = ['v40', 'v29', 'lapa', 'gtex',
                   'encode_cage', 'fantom_cage',
                   'encode_rampage',
                   'pls', 'pels', 'dels', 'ca_h3k4me3',
                   'lrgasp_cage', 'encode_procap']
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
                 input.lrgasp_cage,
                 input.encode_procap]
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
                      end_mode='tes')[0],
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
        gtex = config['gtex']['cerberus']['ics']
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
                    ccre_type='ca_h3k4me3')[0]
    resources:
        mem_gb = 1,
        threads = 1
    params:
        add_ends = [True, True, True,
                  False, False, False, False, False],
        refs = [True, True, False,
                  False, False, False, False, False],
        sources = ['vM25', 'vM21', 'lapa',
                 'fantom_cage', 'pls', 'pels', 'dels', 'h3k4me3']
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
               input.dels]
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

rule all_cerberus:
    input:
        expand(config['lr']['cerberus']['agg_ends'],
               species=species,
               end_mode=end_modes),
        expand(config['lr']['cerberus']['agg_ics'],
              species=species)
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
