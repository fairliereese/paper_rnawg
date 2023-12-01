use rule cerb_gtf_to_bed as cerb_get_gtf_ends_iq with:
    input:
        gtf = config['lr']['isoquant']['gtf']
    output:
        ends = config['lr']['isoquant']['cerberus']['ends']
    params:
        slack = lambda wc:config['params']['cerberus'][wc.end_mode]['slack'],
        dist = lambda wc:config['params']['cerberus'][wc.end_mode]['dist']

use rule cerb_gtf_to_ics as cerb_get_gtf_ics_iq with:
    input:
        gtf = config['lr']['isoquant']['gtf']
    output:
        ics = config['lr']['isoquant']['cerberus']['ics']

rule cerb_make_ref_iq:
    input:
        tss = lambda wc:expand(config['lr']['isoquant']['cerberus']['ends'],
                               species=wc.species,
                               end_mode='tss')[0],
        tes = lambda wc:expand(config['lr']['isoquant']['cerberus']['ends'],
                              species=wc.species,
                              end_mode='tes')[0],
        ic = lambda wc:expand(config['lr']['isoquant']['cerberus']['ics'],
                              species=wc.species)[0],
        h5 = config['lr']['cerberus']['ca_annot']
    resources:
        threads = 1,
        mem_gb = 64
    params:
        source = 'isoquant_wtc11',
        tss_agg_slack = lambda wc:config['params']['cerberus']['tss']['agg_slack'],
        tes_agg_slack = lambda wc:config['params']['cerberus']['tes']['agg_slack'],
        ref = False,
        add_ends = False
    output:
        h5 = config['lr']['isoquant']['cerberus']['ca']
    run:
        ca = cerberus.read(input.h5)
        source = params.source
        ref = params.ref
        add_ends = params.add_ends
        tss_slack = params.tss_agg_slack
        tes_slack = params.tes_agg_slack

        ca.add_bed(input.tss, add_ends,
                   ref, source,
                   'tss', slack=tss_slack)
        ca.add_bed(input.tes, add_ends,
                   ref, source,
                  'tes', slack=tes_slack)
        ca.add_ics(input.ic,
                  ref, source)
        ca.write(output.h5)

use rule cerb_annot as cerberus_annotate_iq with:
    input:
        gtf = config['lr']['isoquant']['gtf'],
        h5 = config['lr']['isoquant']['cerberus']['ca']
    params:
        source = 'isoquant_wtc11',
        gene_source = lambda wc:config['ref'][wc.species]['new_gtf_ver']
    output:
        h5 = config['lr']['isoquant']['cerberus']['ca_annot']

use rule cerb_gtf_ids as cerb_gtf_ids_iq with:
    input:
        h5 = config['lr']['isoquant']['cerberus']['ca_annot'],
        gtf = config['lr']['isoquant']['gtf']
    params:
        source = 'lapa',
        update_ends = True,
        agg = True
    output:
        gtf = config['lr']['isoquant']['cerberus']['gtf']

rule all_isoquant:
    input:
        expand(config['lr']['isoquant']['cerberus']['gtf'],
               species='human')
        # expand(config['lr']['isoquant']['cerberus']['ics'],
        #        species='human'),
        # expand(config['lr']['isoquant']['cerberus']['ends'],
        #       species='human',
        #       end_mode=['tes', 'tss']),
