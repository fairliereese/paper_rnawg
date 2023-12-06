# wildcard_constraints:
#     species='|'.join([re.escape(x) for x in lr_df.species.unique().tolist()]),


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


def format_isoquant_ab(ab, lib_meta, ofile):
    """
    Format abundance file from isoquant GFF into a format
    that Cerberus can deal with
    """
    meta = pd.read_csv(lib_meta, sep='\t')
    df = pd.read_csv(ab, sep='\t')
    df.head()
    df['annot_transcript_id'] = df['#feature_id']
    df['annot_transcript_name'] = df['annot_transcript_id']
    df['transcript_ID'] = df['annot_transcript_id']
    df.columns = [c.split('.')[0] for c in df.columns]
    df.drop('#feature_id', axis=1, inplace=True)
    m = dict([(entry['ENCODE_alignments_id'], entry['dataset']) for ind, entry in meta[['dataset', 'ENCODE_alignments_id']].iterrows()])
    df.columns = [m[c] if 'ENCFF' in c else c for c in df.columns]

    df.to_csv(ofile, sep='\t', index=False)

rule fmt_ab_iq:
    input:
        ab = config['lr']['isoquant']['novel_ab'],
        lib_meta = config['lr']['meta']
    resources:
        threads = 1,
        mem_gb = 8
    output:
        ab_fmt = config['lr']['isoquant']['ab_fmt'],
    run:
        format_isoquant_ab(input.ab, input.lib_meta, output.ab_fmt)

use rule cerb_ab_ids as cerb_ab_ids_iq with:
    input:
        h5 = config['lr']['isoquant']['cerberus']['ca_annot'],
        ab = config['lr']['isoquant']['ab_fmt']
    params:
        source = 'isoquant_wtc11',
        agg = True
    output:
        ab = config['lr']['isoquant']['cerberus']['ab']

###### Tracks
rule tracks_gene_pred_iq:
    input:
        ifile = config['lr']['isoquant']['cerberus']['gtf']
    resources:
        mem_gb = 16,
        threads = 1
    output:
        ofile = temporary(config['lr']['isoquant']['tracks']['gp'])
    shell:
        """
        gtfToGenePred -genePredExt {input.ifile} {output.ofile}
        """

rule tracks_big_gene_pred_iq:
    input:
        ifile = config['lr']['isoquant']['tracks']['gp']
    resources:
        mem_gb = 16,
        threads = 1
    output:
        ofile = temporary(config['lr']['isoquant']['tracks']['bgp'])
    shell:
        """
        genePredToBigGenePred {input.ifile} {output.ofile}
        """

use rule bed_sort as tracks_bed_sort_iq with:
    input:
        ifile = config['lr']['isoquant']['tracks']['bgp']
    output:
        ofile = temporary(config['lr']['isoquant']['tracks']['bgp_sort'])

rule tracks_filt_iq:
    input:
        ifile = config['lr']['isoquant']['tracks']['bgp']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        ofile = config['lr']['isoquant']['tracks']['bgp_sort_filt']
    shell:
        """
        grep -v chrM {input.ifile} > {output.ofile}
        """

rule tracks_bigbed_iq:
    input:
        ifile = config['lr']['isoquant']['tracks']['bgp_sort_filt'],
        # as_file = config['ref']['tracks']['as'],
        as_file = config['ref']['ucsc']['as'],
        chrom_sizes = config['ref']['talon']['chrom_sizes']
    resources:
        mem_gb = 16,
        threads = 1
    output:
        ofile = config['lr']['isoquant']['tracks']['bb']
    shell:
        """
        bedToBigBed -type=bed12+8 -tab -as={input.as_file} {input.ifile} {input.chrom_sizes} {output.ofile}
        """

rule all_isoquant:
    input:
        expand(config['lr']['isoquant']['tracks']['bb'],
               species='human'),
        expand(config['lr']['isoquant']['cerberus']['gtf'],
               species='human'),
        expand(config['lr']['isoquant']['cerberus']['ab'],
              species='human')
        # expand(config['lr']['isoquant']['cerberus']['ics'],
        #        species='human'),
        # expand(config['lr']['isoquant']['cerberus']['ends'],
        #       species='human',
        #       end_mode=['tes', 'tss']),
