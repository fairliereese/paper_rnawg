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
        trips = config['lr']['cerberus']['ca_triplets'],
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

def get_fusion_sample_t_coords(ab, gtf, min_tpm, sample, species, ofile):
    """
    Get genomic start / stop coords for transcripts expressed in
    a given sample
    """

    df = pd.read_csv(ab, sep='\t')
    tids = df.loc[df.gene_novelty=='Fusion', 'annot_transcript_id'].tolist()
    print(len(tids))
    df = get_det_table(df,
                     groupby='sample',
                     how='iso',
                     min_tpm=min_tpm,
                     species=species)
    df = df.transpose()
    tids2 = df.loc[df[sample]==True].index.tolist()
    tids = list(set(tids)&set(tids2))

    gtf_df = pr.read_gtf(gtf).df
    gtf_df = gtf_df.loc[gtf_df.transcript_id.isin(tids)]

    gtf_df = gtf_df.loc[gtf_df.Feature=='transcript']
    gtf_df = gtf_df[['Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'transcript_id',]]

    gtf_df.to_csv(ofile, sep='\t', index=False)

rule get_fusion_sample_t_coords:
    input:
        filt_ab = config['lr']['cerberus']['filt_ab'],
        gtf = config['lr']['cerberus']['gtf']
    resources:
        mem_gb = 16,
        threads = 1
    params:
        min_tpm = 1
    output:
        tsv = config['lr']['analysis']['fusion_t_coords']
    run:
        get_fusion_sample_t_coords(input.filt_ab,
            input.gtf,
            params.min_tpm,
            wildcards.sample,
            wildcards.species,
            output.tsv)

rule all_analysis:
    input:
        list(set(expand(expand(config['lr']['analysis']['fusion_t_coords'],
               zip,
               species=lr_df['species'].tolist(),
               sample=lr_df['sample'].tolist(),
               allow_missing=True),
               obs_col='sample'))),
        expand(config['lr']['analysis']['major_isos'],
               species=species,
               obs_col='sample'),
        expand(config['lr']['cerberus']['ca_triplets'],
               species=species,
               obs_col='sample')
