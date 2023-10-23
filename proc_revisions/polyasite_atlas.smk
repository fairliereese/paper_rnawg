import pyranges as pr

species = ['human', 'mouse']

use rule dl as dl_polyasite_atlas with:
    params:
        link = lambda wc:config['polya_atlas'][wc.species]['link']
    output:
        out = temporary(config['polya_atlas']['bed_gz'])

use rule gunzip as gz_polyasite_atlas with:
    input:
        gz = config['polya_atlas']['bed_gz']
    output:
        out = config['polya_atlas']['bed']

rule format_polya_atlas:
    input:
        bed = config['polya_atlas']['bed']
    resources:
        mem_gb = 2,
        threads = 1
    output:
        fmt_bed = config['polya_atlas']['bed_formatted']
    run:
        # remove noncanonical chromosomes, add 'chr' prefix
        df = pr.read_bed(input.bed, as_df=True)
        df.Chromosome = df.Chromosome.astype(str)
        drop_chr_prefs = ['GL', 'KI']
        for pref in drop_chr_prefs:
            df = df.loc[~df.Chromosome.str.startswith(pref)]
        df.Chromosome = 'chr'+df.Chromosome
        df = pr.PyRanges(df)
        df.to_bed(output.fmt_bed)

rule all_polyasite:
    input:
        expand(config['polya_atlas']['bed_formatted'],
               species=species)
