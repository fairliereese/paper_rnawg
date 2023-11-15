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


rule calc_triplets:
    input:
        swan_file = config['lr']['swan']['sg'],
        h5 = config['lr']['cerberus']['ca_annot'],
        filt_ab = config['lr']['cerberus']['filt_ab'],
        major_isos = config['lr']['analysis']['major_isos'],
    params:
        min_tpm = 1,
        gene_subset = 'polya',
        obs_col = 'sample'
    resources:
        threads = 1,
        mem_gb = 64
    output:
        trips = config['lr']['cerberus']['ca_triplets']
        tsv = config['lr']['analysis']['triplets']
    run:
        if wildcards.species == 'human':
            calculate_human_triplets(input.swan_file,
                                input.h5,
                                input.filt_ab,
                                input.major_isos,
                                output.trips,
                                output.tsv,
                                obs_col=params.obs_col,
                                min_tpm=params.min_tpm,
                                gene_subset=params.gene_subset)
        elif wildcards.species == 'mouse':
            calculate_mouse_triplets(input.swan_file,
                              input.h5,
                              input.filt_ab,
                              input.major_isos,
                              output.trips,
                              output.tsv,
                              obs_col=params.obs_col,
                              min_tpm=params.min_tpm,
                              gene_subset=params.gene_subset)

rule all_analysis:
    input:
        expand(config['lr']['analysis']['major_isos'],
               species=species,
               obs_col='sample'),
        expand(config['lr']['cerberus']['ca_triplets'],
               species=species,
               obs_col='sample')
