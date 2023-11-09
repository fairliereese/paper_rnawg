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
