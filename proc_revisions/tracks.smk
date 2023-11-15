################################################################################
######################## GTF tracks for each sample ############################
################################################################################
rule tracks_sample_gtf:
    input:
        ab = config['lr']['cerberus']['filt_ab'],
        gtf = config['lr']['cerberus']['gtf']
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
        ifile = config['lr']['tracks']['sample']['bgp']
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
        as_file = config['ref']['ucsc']['as'],
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
        expand(config['lr']['tracks']['sample']['bb'],
               zip,
               species=lr_df['species'].tolist(),
               sample=lr_df['sample'].tolist())
