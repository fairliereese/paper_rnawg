import pyranges as pr
df = pr.read_gtf('/pub/smorabit/references/panTro6/genes/panTro6.ncbiRefSeq.gtf').df
df = df[['Chromosome', 'gene_id', 'gene_name']].drop_duplicates()
df.to_csv('panTro6.chrom_to_gene.tsv', sep='\t', index=False)