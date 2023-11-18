################################################################################
######################## GTF tracks for each sample ############################
################################################################################
rule tracks_sample_gtf:
    input:
        ab = config['lr']['cerberus']['filt_ab'],
        gtf = config['lr']['cerberus']['gtf'],
        ref = config['ref']['cerberus']['new_gtf_g_info']
    params:
        min_tpm = 1
    resources:
        mem_gb = 16,
        threads = 2
    output:
        gtf = config['lr']['tracks']['sample']['gtf']
    run:
        get_sample_gtf(input.ab,
                       input.gtf,
                       params.min_tpm,
                       wildcards.sample,
                       wildcards.species,
                       output.gtf)

rule tracks_gene_pred:
    input:
        ifile = config['lr']['tracks']['sample']['gtf']
    resources:
        mem_gb = 16,
        threads = 1
    output:
        ofile = temporary(config['lr']['tracks']['sample']['gp'])
    shell:
        """
        gtfToGenePred -genePredExt {input.ifile} {output.ofile}
        """

rule tracks_big_gene_pred:
    input:
        ifile = config['lr']['tracks']['sample']['gp']
    resources:
        mem_gb = 16,
        threads = 1
    output:
        ofile = temporary(config['lr']['tracks']['sample']['bgp'])
    shell:
        """
        genePredToBigGenePred {input.ifile} {output.ofile}
        """

rule tracks_add_bgp_data:
    input:
        ifile = rules.tracks_big_gene_pred.output.ofile,
        h5 = config['lr']['cerberus']['ca_triplets'],
        filt_ab = config['lr']['cerberus']['filt_ab'],
        swan_file = config['lr']['swan']['sg'],
        meta_file = config['lr']['meta'],
        ppred_file = [],
        pi_table = config['lr']['mane']['pi_tpm']['triplet'],
        major_isos = config['lr']['analysis']['major_isos']
    resources:
        mem_gb = 16,
        threads = 4
    params:
        min_tpm = 1
    output:
        ofile = config['lr']['tracks']['sample']['bgp_plus']
    run:
        add_bgp_info(input.ifile,
                     wildcards.species,
                     wildcards.obs_col,
                     params.min_tpm,
                     wildcards.sample,
                     input.h5,
                     input.filt_ab,
                     input.swan_file,
                     input.meta_file,
                     input.ppred_file,
                     input.pi_table,
                     input.major_isos,
                     output.ofile)

rule bed_sort:
    resources:
        mem_gb = 64,
        threads = 1
    shell:
        """
        sort -k1,1 -k2,2n {input.ifile} > {output.ofile}
        """

use rule bed_sort as tracks_bed_sort with:
    input:
        ifile = rules.tracks_add_bgp_data.output.ofile
    output:
        ofile = temporary(config['lr']['tracks']['sample']['bgp_sort'])

rule tracks_filt:
    input:
        ifile = config['lr']['tracks']['sample']['bgp_sort']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        ofile = config['lr']['tracks']['sample']['bgp_sort_filt']
    shell:
        """
        grep -v chrM {input.ifile} > {output.ofile}
        """

use rule dl as dl_ucsc_as with:
  params:
    link = config['ref']['ucsc']['as_link']
  output:
    out = config['ref']['ucsc']['as']

rule tracks_bigbed:
    input:
        ifile = config['lr']['tracks']['sample']['bgp_sort_filt'],
        as_file = config['ref']['tracks']['as'],
        # as_file = config['ref']['ucsc']['as'],
        chrom_sizes = config['ref']['talon']['chrom_sizes']
    resources:
        mem_gb = 16,
        threads = 1
    output:
        ofile = config['lr']['tracks']['sample']['bb']
    shell:
        """
        bedToBigBed -type=bed12+8 -tab -as={input.as_file} {input.ifile} {input.chrom_sizes} {output.ofile}
        """
rule all_tracks:
    input:
        list(set(expand(expand(config['lr']['tracks']['sample']['bb'],
               zip,
               species=lr_df['species'].tolist(),
               sample=lr_df['sample'].tolist(),
               allow_missing=True),
               obs_col='sample')))
