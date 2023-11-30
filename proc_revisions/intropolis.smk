################################################################################
############################### Intropolis #####################################
################################################################################
use rule dl as dl_intropolis with:
    params:
        link = config['intropolis']['link']
    output:
        out = temporary(config['intropolis']['tab_gz'])

use rule gunzip as gz_intrpolis with:
    input:
        gz = config['intropolis']['tab_gz']
    output:
        out = temporary(config['intropolis']['tab'])

rule format_intropolis:
    input:
        tab = config['intropolis']['tab']
    resources:
        mem_gb = 32,
        threads = 2
    output:
        bed = config['intropolis']['bed']
    run:
        df = pd.read_csv(input.tab,
                         header=None,
                         usecols=[8,9,10,11],
                         names=['Chromosome', 'Start',
                                'End', 'Strand'],
                         sep='\t')
        df['Start'] = df['Start'] - 1
        df = df.loc[df.Start.notnull()]
        df = pr.PyRanges(df)
        df.to_bed(output.bed)

rule all_intropolis:
    input:
        expand(config['intropolis']['bed'],
               species='human')
