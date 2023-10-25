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
######################### Cerberus aggregation #################################
################################################################################

rule cerb_agg_ends:
    resources:
        threads = 4,
        mem_gb = 64
    run:
        pass
        # TODO

rule cerberus_agg_ics:
    resources:
        mem_gb = 32,
        threads = 2
    run:
        pass
        # TODO

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

#### agg stuff
rule cerb_agg_human_tss_config:
    input:
        v40 = expand(config['ref']['cerberus']['new_ends'],
                     species='human',
                     end_mode='tss'),
        v29 = expand(config['ref']['cerberus']['ends'],
                     species='human',
                     end_mode='tss'),
        lapa = expand(config['lr']['cerberus']['ends'],
                      species='human',
                      end_mode='tss'),
        gtex = expand(config['gtex']['cerberus']['ends'],
                      end_mode='tss'),
        encode_cage = expand(config['cage']['merged'],
                             species='human'),
        fantom_cage = expand(config['fantom']['bed'],
                             species='human'),
        encode_rampage = expand(config['rampage']['merged'],
                                species='human'),
        pls = expand(config['ccre']['bed_format'],
                     species='human',
                     ccre_type='pls'),
        pels = expand(config['ccre']['bed_format'],
                      species='human',
                      ccre_type='pels'),
        dels = expand(config['ccre']['bed_format'],
                      species='human',
                      ccre_type='dels'),
        lrgasp_cage = expand(config['lrgasp_cage']['merged'],
                             species='human'),
        encode_procap = expand(config['procap']['merged'],
                               species='human')
    params:
        add_ends = [True, True, True, True,
                    False, False, False,
                    False, False, False, False, False],
        refs = [True, True, False, False,
                False, False, False,
                False, False, False, False, False],
        sources = ['v40', 'v29', 'lapa', 'gtex',
                   'encode_cage', 'fantom_cage',
                   'encode_rampage',
                   'pls', 'pels', 'dels',
                   'lrgasp_cage', 'encode_procap']
    output:
        cfg = expand(config['lr']['cerberus']['agg_ends_cfg'],
                     species='human',
                     end_mode='tss')
    run:
        files = [input.v40,
                 input.v29,
                 input.lapa_gtf,
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
                     end_mode='tes'),
        v29 = expand(config['ref']['cerberus']['ends'],
                     species='human',
                     end_mode='tes'),
        lapa = expand(config['lr']['cerberus']['ends'],
                      species='human',
                      end_mode='tes'),
        gtex = expand(config['gtex']['cerberus']['ends'],
                      end_mode='tes'),
        pas = expand(config['pas']['ends_formatted'],
                species='human',
                end_mode='tes'),
        atlas = expand(config['polya_atlas']['bed_formatted'],
                       species='human')
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
                     end_mode='tes')
    run:
         files = [input.v40,
                  input.v29,
                  input.lapa_gtf,
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
                     species='human'),
        v29 = expand(config['ref']['cerberus']['ics'],
                     species='human'),
        lapa = expand(config['lr']['cerberus']['ics'],
                     species='human'),
        gtex = config['gtex']['cerberus']['ics']
    params:
        sources = ['v40', 'v29', 'lapa', 'gtex'],
        refs = [True, True, False, False]
    resources:
        mem_gb = 1,
        threads = 1
    output:
        cfg = config['lr']['cerberus']['agg_ics_cfg']
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
                     end_mode='tss'),
        vM25 = expand(config['ref']['cerberus']['new_ends'],
                     species='mouse',
                     end_mode='tss'),
        lapa = expand(config['lr']['cerberus']['ends'],
                    species='mouse',
                    end_mode='tss'),
        fantom_cage = expand(config['fantom']['bed'],
                             species='mouse'),
        pls = expand(config['ccre']['bed_format'],
                     species='mouse',
                     ccre_type='pls'),
        pels = expand(config['ccre']['bed_format'],
                      species='mouse',
                      ccre_type='pels'),
        dels = expand(config['ccre']['bed_format'],
                      species='mouse',
                      ccre_type='dels')
    resources:
        mem_gb = 1,
        threads = 1
    params:
        add_ends = [True, True, True,
                  False, False, False, False],
        refs = [True, True, False,
                  False, False, False, False],
        sources = ['vM25', 'vM21', 'lapa',
                 'fantom_cage', 'pls', 'pels', 'dels']
    output:
        cfg = expand(config['lr']['cerberus']['agg_ends_cfg'],
                   species='mouse',
                   end_mode='tss')
    run:
        files = [input.vM25,
               input.vM21,
               input.lapa_gtf,
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
                     end_mode='tes'),
        vM25 = expand(config['ref']['cerberus']['new_ends'],
                     species='mouse',
                     end_mode='tes'),
        lapa = expand(config['lr']['cerberus']['ends'],
                    species='mouse',
                    end_mode='tes'),
        pas = expand(config['pas']['ends_formatted'],
                species='mouse',
                end_mode='tes'),
        atlas = expand(config['polya_atlas']['bed_formatted'],
                       species='mouse')
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
        cfg = config['cerberus']['tes']['cfg']
    run:
        files = [input.vM25,
                input.vM21,
                input.lapa_gtf,
                input.pas,
                input.atlas]
        df = pd.DataFrame()
        df['fname'] = files
        df['add_ends'] = params.add_ends
        df['refs'] = params.refs
        df['sources'] = params.sources
        df.to_csv(output.cfg, sep=',', header=None, index=False)


rule cerb_agg_mouse_ic_config:
    input:
        vM25 = expand(config['ref']['cerberus']['new_ics'],
                     species='mouse'),
        vM21 = expand(config['ref']['cerberus']['ics'],
                     species='mouse'),
        lapa = expand(config['lr']['cerberus']['ics'],
                     species='mouse')
    params:
        sources = ['vM25', 'vM21', 'lapa'],
        refs = [True, True, False]
    resources:
        mem_gb = 1,
        threads = 1
    output:
        cfg = config['lr']['cerberus']['agg_ics_cfg']
    run:
        files = [input.vM25, input.vM21, input.lapa]
        refs = params.refs
        sources = params.sources
        df = pd.DataFrame()
        df['fname'] = files
        df['refs'] = refs
        df['source'] = sources
        df.to_csv(output.cfg, header=None, index=False, sep=',')


rule all_cerberus:
    input:
        expand(config['ref']['cerberus']['ends'],
              species=species,
              end_mode=end_modes),
        expand(config['ref']['cerberus']['new_ends'],
            species=species,
            end_mode=end_modes),
        expand(config['ref']['cerberus']['ics'],
              species=species),
        expand(config['ref']['cerberus']['new_ics'],
            species=species),
        expand(config['ref']['cerberus']['ends'],
              species=species,
              end_mode=end_modes),
        expand(config['lr']['cerberus']['ics'],
            species=species)
