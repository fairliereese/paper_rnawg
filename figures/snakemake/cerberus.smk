import cerberus

################################################################################
############################## Cerberus ########################################
################################################################################
rule cerb_gtf_to_bed:
    resources:
        mem_gb = 64,
        threads = 1
    run:
        cerberus.gtf_to_bed(input.gtf,
                            wildcards.end_mode,
                            output.ends,
                            dist=params.dist,
                            slack=params.slack)

rule cerb_gtf_to_ics:
    resources:
        mem_gb = 64,
        threads = 1
    run:
        cerberus.gtf_to_ics(input.gtf,
                            output.ics)

rule cerb_agg_ends:
    resources:
      threads = 4,
      mem_gb = 32
    run:
        refs = [params.refs for i in range(len(input.files))]
        add_ends = [params.add_ends for i in range(len(input.files))]
        cerberus.agg_ends(input.files,
                          add_ends,
                          refs,
                          params.sources,
                          wildcards.end_mode,
                          params.slack,
                          output.ends)

rule cerb_agg_ics:
  resources:
    threads = 4,
    mem_gb = 32
  run:
      refs = [params.refs for i in range(len(input.files))]
      cerberus.agg_ics(input.files,
                        refs,
                        params.sources,
                        output.ics)

rule cerb_write_ref:
    resources:
        threads = 4,
        mem_gb = 64
    run:
        cerberus.write_reference(input.tss,
                                 input.tes,
                                 input.ic,
                                 output.h5)


################################################################################
######################### Cerberus annot + ID replacement ######################
################################################################################

rule cerb_annot:
    resources:
        mem_gb = 64,
        threads = 16
    run:
        cerberus.annotate_transcriptome(input.gtf,
                                        input.h5,
                                        params.source,
                                        params.gene_source,
                                        output.h5)

rule cerb_gtf_ids:
    resources:
        mem_gb = 64,
        threads = 16
    run:
        cerberus.replace_gtf_ids(input.h5,
                                 input.gtf,
                                 params.source,
                                 True,
                                 True,
                                 output.gtf)

rule cerb_ab_ids:
    resources:
        mem_gb = 64,
        threads = 16
    run:
        cerberus.replace_ab_ids(input.ab,
                                input.h5,
                                params.source,
                                True,
                                output.ab)
