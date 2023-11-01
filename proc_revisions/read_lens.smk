# lr meta
meta_df = get_meta_df(config, species)

################################################################################
############################## Diane's stuff ###################################
################################################################################
rule dl_lr_bam:
    resources:
        mem_gb = 32,
        threads = 8
    params:
        encid = lambda w:get_encid_from_dataset(w.dataset,
                                                meta_df,
                                                'bam')
    output:
        bam = temporary(config['lr']['bam'])
    shell:
        "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.bam -O {output.bam}"

rule dl_lr_label_bam:
    resources:
        mem_gb = 32,
        threads = 8
    params:
        encid = lambda w:get_encid_from_dataset(w.dataset,
                                                meta_df,
                                                'label_bam')
    output:
        bam = temporary(config['lr']['label_bam'])
    shell:
        "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.bam -O {output.bam}"


rule dl_lr_fastq:
    resources:
        mem_gb = 32,
        threads = 8
    params:
        encid = lambda w:get_encid_from_dataset(w.dataset,
                                                meta_df,
                                                'fastq')
    output:
        fastq = temporary(config['lr']['fastq_gz'])
    shell:
        "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.fastq.gz -O {output.fastq}"

rule get_lr_read_lens:
    input:
        bams = expand(config['data']['bam'],
                      species='human',
                      dataset=lambda w:get_col_from_meta_df(wc, col='dataset')),
        fastqs = expand(config['data']['fastq_gz'],
                        species='human',
                        dataset=lambda w:get_col_from_meta_df(wc, col='dataset'))
    resources:
        mem_gb = 32,
        threads = 8
    output:
        tsv = config['data']['read_len_meta']
    run:
        get_lr_read_lens(input.bams, input.fastqs, output.tsv)
