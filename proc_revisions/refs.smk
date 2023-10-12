include: 'download.smk'
configfile: 'config.yml'

################################################################################
######################## Ref. download / proc ##################################
################################################################################
use rule dl as dl_fa with:
    params:
        link = lambda wc:config['ref'][wc.species]['fa_link']
    output:
        out = temporary(config['ref']['fa_gz'])

use rule gunzip as gz_fa with:
    input:
        gz = config['ref']['fa_gz']
    output:
        out = config['ref']['fa']

use rule dl as dl_annot with:
    params:
        link = lambda wc:config['ref'][wc.species]['gtf_link']
    output:
        out = temporary(config['ref']['gtf_gz'])

use rule gunzip as gz_annot with:
    input:
        gz = config['ref']['gtf_gz']
    output:
        out = config['ref']['gtf']

use rule dl as dl_new_annot with:
    params:
        link = lambda wc:config['ref'][wc.species]['new_gtf_link']
    output:
        out = temporary(config['ref']['new_gtf_gz'])

use rule gunzip as gz_new_annot with:
    input:
        gz = config['ref']['new_gtf_gz']
    output:
        out = config['ref']['new_gtf']

#################################################################################
####################### Ref. w/ spikes for TALON ################################
################################################################################

rule mkref_gtf_with_spike:
    input:
        sirv = config['ref']['spike']['sirv_gtf'],
        ercc = config['ref']['spike']['ercc_gtf'],
        gtf = config['ref']['gtf']
    resources:
        mem_gb = 16,
        threads = 2
    output:
        all = config['ref']['talon']['gtf']
    shell:
        """
        cat {input.gtf} > {output.all}
        cat {input.sirv} >> {output.all}
        cat {input.ercc} >> {output.all}
        """

rule mkref_spike_fa:
    input:
        sirv = config['ref']['spike']['sirv_fa'],
        ercc = config['ref']['spike']['ercc_fa'],
        fa = config['ref']['fa']
    resources:
        threads = 1,
        mem_gb = 4
    output:
        all = config['ref']['talon']['fa']
    shell:
        """
        cat {input.fa} >> {output.all}
        cat {input.sirv} >> {output.all}
        cat {input.ercc} >> {output.all}
        """

rule get_chrom_sizes:
    input:
        fa = config['ref']['talon']['fa']
    resources:
        threads = 1,
        mem_gb = 8
    output:
        chrom_sizes = config['ref']['talon']['chrom_sizes']
    shell:
        """
        faidx {input.fa} -i chromsizes > {output.chrom_sizes}
        """

rule get_utr_fix_gtf:
    input:
        gtf = config['ref']['talon']['gtf']
    resources:
        threads = 1,
        mem_gb = 8
    output:
        gtf = config['ref']['lapa']['gtf_utr']
    shell:
        """
        gencode_utr_fix \
            --input_gtf {input.gtf} \
            --output_gtf {output.gtf}
        """
