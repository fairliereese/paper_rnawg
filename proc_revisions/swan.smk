rule swan_init:
    input:
        annot = config, # todo
        ab = config['lr']['cerberus']['filt_ab'],
        gene_ab = config['lr']['talon']['ab'],
        gtf = config['lr']['cerberus']['gtf'],
        meta = config['lr']['meta']
    params:
        prefix = config['swan']['sg'].replace('.p', '')
    resources:
        mem_gb = 64,
        threads = 1
    output:
        sg = config['swan']['sg']
    run:
        make_sg(input, params, wildcards)
