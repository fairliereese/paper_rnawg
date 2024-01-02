use rule dl as dl_gtex with:
    params:
        link = config['gtex']['gtf_link']
    output:
        out = temporary(config['gtex']['gtf_gz'])

use rule gunzip as gz_gtex with:
    input:
        gz = config['gtex']['gtf_gz']
    output:
        out = config['gtex']['gtf']

# create gtex gtf only with known genes
rule filter_gtex_gtf:
    input:
        gtf = config['gtex']['gtf']
    resources:
        mem_gb = 8,
        threads = 1
    output:
        filt_gtf = config['gtex']['filt_gtf']
    run:
        filter_gtex_gtf(input.gtf, output.filt_gtf)

use rule cerb_gtf_to_bed as cerb_get_gtf_ends with:
    input:
        gtf = config['gtex']['filt_gtf']
    output:
        ends = config['gtex']['cerberus']['ends']
    params:
        slack = lambda wc:config['params']['cerberus'][wc.end_mode]['slack'],
        dist = lambda wc:config['params']['cerberus'][wc.end_mode]['dist']

use rule cerb_gtf_to_ics as cerb_get_gtf_ics with:
    input:
        gtf = config['gtex']['filt_gtf']
    output:
        ics = config['gtex']['cerberus']['ics']

use rule cerb_gtf_ids as cerb_gtf_ids_new_gtex with:
    input:
        h5 = config['lr']['cerberus']['ca_annot'],
        gtf = config['gtex']['filt_gtf']
    params:
        source = lambda wc:config['ref'][wc.species]['new_gtf_ver'],
        update_ends = True,
        agg = True
    output:
        gtf = config['gtex']['cerberus']['gtf']

use rule dl as dl_gtex_ab with:
    params:
        link = config['gtex']['ab_link']
    output:
        out = temporary(config['gtex']['ab_gz'])

use rule gunzip as gz_gtex_ab with:
    input:
        gz = config['gtex']['ab_gz']
    output:
        out = config['gtex']['ab']

def format_gtex_abundance(ifile, ofile):
    df = pd.read_csv(ifile, sep='\t')
    df['annot_transcript_id'] = df.transcript
    df['annot_transcript_name'] = df.transcript
    df['transcript_ID'] = df.transcript_ID
    df.drop('transcript', axis=1, inplace=True)
    df.to_csv(ofile, sep='\t', index=False)

rule gtex_fmt_ab:
    input:
        ab = config['gtex']['ab']
    resources:
        threads = 1,
        mem_gb = 4
    output:
        ab = config['gtex']['ab_fmt']
    run:
        format_gtex_abundance(input.ab, output.ab)

use rule cerb_ab_ids as cerb_ab_ids_gtex with:
    input:
        h5 = config['lr']['cerberus']['ca_annot'],
        ab = config['gtex']['ab_fmt']
    params:
        source = 'gtex',
        agg = True
    output:
        ab = config['gtex']['cerberus']['ab']

rule all_gtex:
    input:
        expand(config['gtex']['cerberus']['gtf'],
               species='human'),
        expand(config['gtex']['cerberus']['ab'], species='human')
