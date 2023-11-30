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

# def get_tss_ca_add_settings(wc, how):
#     """
#     Get files or sources for additional tss support
#
#     Parameters:
#         how (str): {'file', 'source'}
#     """
#
#     files = []
#     sources = []
#
#     # procap bidirectional
#     files += expand(config['procap']['bed'],
#                     procap_dataset=procap_meta.dataset.unique().tolist(),
#                     species=wc.species,
#                     output_type='bidirectional_peaks')
#     sources += [d+'_bi_procap' for d in procap_meta.dataset.unique().tolist()]
#
#     # procap unidirectional
#     files += expand(config['procap']['format_uni_bed'],
#                     procap_dataset=procap_meta.dataset.unique().tolist(),
#                     species=wc.species,
#                     output_type='unidirectional_peaks')
#     sources += [d+'_uni_procap' for d in procap_meta.dataset.unique().tolist()]
#
#     # lrgasp cage
#     files += expand(config['lrgasp_cage']['bed'],
#                     lrgasp_cage_dataset=lrgasp_cage_meta.dataset.unique().tolist(),
#                     species=wc.species)
#     sources += [d+'_lrgasp_cage' for d in lrgasp_cage_meta.dataset.unique().tolist()]
#
#     if how == 'file':
#         return files
#     elif how == 'source':
#         return sources

rule cerb_add_isoquant_ends:
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

        ca.add_bed(bed, add_ends,
                   ref, source,
                   'tss', slack=tss_slack)
        ca.add_bed(bed, add_ends,
                   ref, source,
                  'tes', slack=tes_slack)


        ca.write(output.h5)

rule all_isoquant:
    input:
        expand(config['lr']['isoquant']['cerberus']['ics'],
               species='human'),
        expand(config['lr']['isoquant']['cerberus']['ends'],
              species='human',
              end_mode=['tes', 'tss']),
