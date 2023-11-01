rule get_transcript_novelties:
    input:
        annot = rules.cerberus_annotate_lr.output.h5,
        filt_ab = rules.cerb_ab_ids_lr.output.ab,
        t_meta = rules.t_info_lr.output.o
    params:
        min_tpm = 1,
        gene_subset = 'polya',
        ver = 'v40_cerberus'
    resources:
        mem_gb = 8,
        threads = 1
    output:
        novelties = config['lr']['suppa']['novelties']
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
        nov = config['lr']['suppa']['novelties'],
        filt_ab = rules.cerb_ab_ids_lr.output.ab
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
#
# rule suppa_generate_events:
#     resources:
#         mem_gb = 16,
#         threads = 1
#     shell:
#         """
#         suppa.py generateEvents \
#             -i {input.gtf} \
#             -o {params.opref} \
#             -f ioe \
#             -e SE SS MX RI FL
#         """
#
# use rule suppa_generate_events as ge_cerb with:
#     input:
#         gtf = config['data']['suppa']['gtf']
#     params:
#         opref = config['data']['suppa']['events']['A3'].rsplit('_', maxsplit=2)[0]
#         # opref = config['data']['suppa']['events']
#     output:
#         config['data']['suppa']['events']['A3'],
#         config['data']['suppa']['events']['A5'],
#         config['data']['suppa']['events']['AF'],
#         config['data']['suppa']['events']['AL'],
#         config['data']['suppa']['events']['MX'],
#         config['data']['suppa']['events']['RI'],
#         config['data']['suppa']['events']['SE'],
#
# use rule suppa_generate_events as ge_gtex with:
#     input:
#         gtf = config['ref']['gtex_gtf']
#     params:
#         opref = config['data']['suppa']['gtex']['events']['A3'].rsplit('_', maxsplit=2)[0]
#     output:
#         config['data']['suppa']['gtex']['events']['A3'],
#         config['data']['suppa']['gtex']['events']['A5'],
#         config['data']['suppa']['gtex']['events']['AF'],
#         config['data']['suppa']['gtex']['events']['AL'],
#         config['data']['suppa']['gtex']['events']['MX'],
#         config['data']['suppa']['gtex']['events']['RI'],
#         config['data']['suppa']['gtex']['events']['SE']
#
# rule suppa_psi:
#     resources:
#         mem_gb = 16,
#         threads = 1
#     shell:
#         """
#         suppa.py psiPerEvent \
#             --ioe-file {input.ioe} \
#             --expression-file {input.filt_ab} \
#             --o {params.opref}
#         """
#
# use rule suppa_psi as psi_cerb with:
#     input:
#         ioe = lambda w:config['data']['suppa']['events'][w.event],
#         filt_ab = config['data']['suppa']['filt_ab']
#     params:
#         opref = config['data']['suppa']['psi'].rsplit('.psi', maxsplit=1)[0]
#     output:
#         out = config['data']['suppa']['psi']
#
# rule get_gtex_cerb_ids:
#     input:
#         cerb_annot = config['data']['cerb_annot']
#     resources:
#         mem_gb = 2,
#         threads = 1
#     output:
#         cerb_ids = config['data']['suppa']['gtex']['cerb_ids']
#     run:
#         get_gtex_cerberus_ids(input.cerb_annot, output.cerb_ids)
#
#
# rule get_cerberus_psi:
#     input:
#         filt_ab = config['data']['filt_ab']
#     params:
#         min_tpm = 1,
#         gene_subset = 'polya'
#     resources:
#         threads = 1,
#         mem_gb = 64
#     output:
#         ofile = config['data']['psi']
#     run:
#         get_cerberus_psi(input.filt_ab,
#                          params.min_tpm,
#                          params.gene_subset,
#                          output.ofile)
#
# e_map = {'tss': 'AF',
#       'tes': 'AL',
#       'AL': 'tes',
#       'AF': 'tss'}
# rule get_matching_events:
#     input:
#         cerb_file = config['data']['psi'],
#         suppa_file = config['data']['suppa']['psi'],
#         lib_meta = config['data']['meta']
#     params:
#         kind = lambda w:e_map[w.event]
#     output:
#         ofile = config['data']['suppa']['matching_events']
#     run:
#         get_cerb_suppa_matching_events(input.cerb_file,
#                                        input.suppa_file,
#                                        output.ofile,
#                                        input.lib_meta,
#                                        params.kind)
