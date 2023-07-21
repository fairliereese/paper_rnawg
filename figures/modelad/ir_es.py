import swan_vis as swan

sg = swan.read('swan_modelad.p')
#es_df = sg.find_es_genes(verbose=True)
#es_df.to_csv('novel_es_events.tsv', sep='\t')
ir_df = sg.find_ir_genes(verbose=True)
ir_df.to_csv('novel_ir_events.tsv', sep='\t')
