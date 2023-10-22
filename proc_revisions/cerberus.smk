import cerberus

end_modes = ['tss', 'tes']

################################################################################
#################### Get triplet features from GTF #############################
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

################################################################################
######################### Cerberus aggregation #################################
################################################################################

rule cerb_agg_ends:
    resources:
        threads = 4,
        mem_gb = 64
    run:
        pass
        # TODO

################################################################################
######################### Cerberus annotation ##################################
################################################################################

rule cerberus_agg_ics:
    resources:
        mem_gb = 32,
        threads = 2
    run:
        # only aggregating 2 things at a time
        refs = [params.ref for i in range(2)]
        cerberus.agg_ics([input.ref_ics, input.ics],
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
                                params.update_ends,
                                params.agg,
                                output.gtf)

rule cerb_ab_ids:
    resources:
        mem_gb = 64,
        threads = 16
    run:
        cerberus.replace_ab_ids(input.ab,
                                input.h5,
                                params.source,
                                params.agg,
                                output.ab)
