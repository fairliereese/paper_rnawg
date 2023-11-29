rule get_transcript_novelties:
    input:
        annot = rules.cerberus_annotate_lr.output.h5,
        filt_ab = config['lr']['cerberus']['filt_ab'],
        t_meta = rules.t_info_lr.output.o
    params:
        min_tpm = 1,
        gene_subset = 'polya',
        ver = 'v40_cerberus'
    resources:
        mem_gb = 8,
        threads = 1
    output:
        novelties = config['lr']['cerberus']['novelties']
    run:
        get_transcript_novelties(input.annot,
                                 input.filt_ab,
                                 input.t_meta,
                                 params.min_tpm,
                                 params.gene_subset,
                                 params.ver,
                                 output.novelties)

rule preproc_suppa:
    input:
        gtf = rules.cerb_gtf_ids_lr.output.gtf,
        nov = config['lr']['cerberus']['novelties'],
        filt_ab = config['lr']['cerberus']['filt_ab']
    params:
        temp_filt_ab_nov = lambda wc:expand('temp_{species}_filt_ab_nov.tsv',
                                        species=wc.species)[0],
        temp_filt_ab_nov_2 = lambda wc:expand('temp_{species}_filt_ab_nov_2.tsv',
                                        species=wc.species)[0]

    output:
        filt_ab = config['lr']['suppa']['filt_ab'],
        gtf = config['lr']['suppa']['gtf']
    conda:
        'seurat'
    shell:
        """
        sed 's/,/_/g' {input.gtf} > {output.gtf}
        Rscript ../scripts/filter_abundance_mat_by_novelty.R \
            --novelties {input.nov} \
            --filtab {input.filt_ab} \
            --ofile {params.temp_filt_ab_nov}
        cut -f4,12- {params.temp_filt_ab_nov} | sed 's/,/_/g' > {params.temp_filt_ab_nov_2}
        python ../scripts/suppa_format_rm_colname.py {params.temp_filt_ab_nov_2}
        tail -c +2 {params.temp_filt_ab_nov_2} > {output.filt_ab}
        """

rule suppa_generate_events:
    resources:
        mem_gb = 16,
        threads = 1
    shell:
        """
        suppa.py generateEvents \
            -i {input.gtf} \
            -o {params.opref} \
            -f ioe \
            -e SE SS MX RI FL
        """

use rule suppa_generate_events as ge_cerb with:
    input:
        gtf = config['lr']['suppa']['gtf']
    params:
        opref = config['lr']['suppa']['events']['A3'].rsplit('_', maxsplit=2)[0]
        # opref = config['lr']['suppa']['events']
    output:
        config['lr']['suppa']['events']['A3'],
        config['lr']['suppa']['events']['A5'],
        config['lr']['suppa']['events']['AF'],
        config['lr']['suppa']['events']['AL'],
        config['lr']['suppa']['events']['MX'],
        config['lr']['suppa']['events']['RI'],
        config['lr']['suppa']['events']['SE']

use rule suppa_generate_events as ge_gtex with:
    input:
        gtf = config['gtex']['gtf']
    params:
        opref = config['gtex']['suppa']['events']['A3'].rsplit('_', maxsplit=2)[0]
    output:
        config['gtex']['suppa']['events']['A3'],
        config['gtex']['suppa']['events']['A5'],
        config['gtex']['suppa']['events']['AF'],
        config['gtex']['suppa']['events']['AL'],
        config['gtex']['suppa']['events']['MX'],
        config['gtex']['suppa']['events']['RI'],
        config['gtex']['suppa']['events']['SE']

rule suppa_psi:
    resources:
        mem_gb = 16,
        threads = 1
    shell:
        """
        suppa.py psiPerEvent \
            --ioe-file {input.ioe} \
            --expression-file {input.filt_ab} \
            --o {params.opref}
        """

use rule suppa_psi as psi_cerb with:
    input:
        ioe = lambda w:config['lr']['suppa']['events'][w.event],
        filt_ab = config['lr']['suppa']['filt_ab']
    params:
        opref = config['lr']['suppa']['psi'].rsplit('.psi', maxsplit=1)[0]
    output:
        out = config['lr']['suppa']['psi']

rule get_gtex_cerb_ids:
    input:
        cerb_annot = config['lr']['cerberus']['ca_triplets']
    resources:
        mem_gb = 2,
        threads = 1
    output:
        cerb_ids = config['gtex']['suppa']['cerb_ids']
    run:
        get_gtex_cerberus_ids(input.cerb_annot, output.cerb_ids)


rule get_cerberus_psi:
    input:
        filt_ab = config['lr']['cerberus']['filt_ab']
    params:
        min_tpm = 1,
        gene_subset = 'polya'
    resources:
        threads = 1,
        mem_gb = 64
    output:
        ofile = config['lr']['cerberus']['psi']
    run:
        get_cerberus_psi(input.filt_ab,
                         params.min_tpm,
                         params.gene_subset,
                         output.ofile)

e_map = {'tss': 'AF',
      'tes': 'AL',
      'AL': 'tes',
      'AF': 'tss'}
events_2 = [item for key, item in e_map.dict()]

rule get_matching_events:
    input:
        cerb_file = config['lr']['cerberus']['psi'],
        suppa_file = config['lr']['suppa']['psi'],
        lib_meta = config['lr']['meta']
    params:
        kind = lambda w:e_map[w.event]
    output:
        ofile = config['lr']['suppa']['matching_events']
    run:
        get_cerb_suppa_matching_events(input.cerb_file,
                                       input.suppa_file,
                                       output.ofile,
                                       input.lib_meta,
                                       params.kind)

rule all_suppa:
    input:
        expand(config['lr']['suppa']['matching_events'],
               species='human',
               events=events_2)
