import swan_vis as swan

def make_sg(input, params, wildcards):

    # initialize
    sg = swan.SwanGraph()
    sg.add_annotation(input.annot)
    sg.add_transcriptome(input.gtf, include_isms=True)
    sg.save_graph(params.prefix)

    sg.add_abundance(input.ab)
    # sg.add_abundance(input.gene_ab, how='gene')
    sg.save_graph(params.prefix)

    # # add metadata and add colors
    # sg.add_metadata(input.meta)
    # c_dict, order = get_biosample_colors(wildcards.species)
    # sg.set_metadata_colors('sample', c_dict)
    #
    # # human only settings
    # if wildcards.species == 'human':
    #     c_dict, order = get_ad_colors()
    #     sg.set_metadata_colors('health_status', c_dict)
    # # save
    # sg.save_graph(params.prefix)

rule swan_gene_ab_add_stable_gid:
    input:
        ab = rules.talon_ab_full.output.ab
    resources:
        mem_gb = 24,
        threads = 1
    output:
        ab = config['lr']['talon']['ab_stable_gid']
    run:
        df = pd.read_csv(input.ab, sep='\t')
        df['gid_stable'] = cerberus.get_stable_gid(df, 'annot_gene_id')
        df['annot_gene_id'] = df['gid_stable']
        df.drop('gid_stable', axis=1, inplace=True)
        df.to_csv(output.ab, sep='\t', index=False)

rule swan_init:
    input:
        annot = config['ref']['cerberus']['new_gtf'],
        ab = rules.cerb_ab_ids_lr.output.ab,
        gene_ab = rules.swan_gene_ab_add_stable_gid.output.ab,
        sgtf = rules.cerb_gtf_ids_lr.output.gtf,
        meta = config['lr']['meta']
    params:
        prefix = config['lr']['swan']['sg'].replace('.p', '')
    resources:
        mem_gb = 64,
        threads = 1
    output:
        sg = config['lr']['swan']['sg']
    run:
        make_sg(input, params, wildcards)

rule all_swan:
    input:
        expand(rules.swan_init.output, species=species)
