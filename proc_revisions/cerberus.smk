import cerberus

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
