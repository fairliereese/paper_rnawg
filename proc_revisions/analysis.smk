from utils import *

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
        ofile = config['lr']['analysis']['major_isos']
    run:
        get_major_isos(input.sg,
                       input.filt_ab,
                       wildcards.obs_col,
                       wildcards.species,
                       output.ofile,
                       min_tpm=params.min_tpm,
                       gene_subset=params.gene_subset)

rule all_analysis:
    input:
        expand(config['lr']['analysis']['major_isos'],
               species=species,
               obs_col='sample')
