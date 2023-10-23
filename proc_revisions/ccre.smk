# ccre stuff
ccre_df = pd.read_csv('ccre_config.tsv', sep='\t')
ccre_types = ccre_df.ccre_type.unique().tolist()
species = ['human', 'mouse']

wildcard_constraints:
    ccre_type='|'.join([re.escape(x) for x in ccre_df.ccre_type.unique().tolist()]),

def get_ccre_link(wc, ccre_df):
    temp = ccre_df.copy(deep=True)
    return ccre_df.loc[(ccre_df.species==wc.species)&\
                       (ccre_df.ccre_type==wc.ccre_type), 'link'].values[0]

use rule dl as dl_ccre with:
    params:
        link = lambda wc:get_ccre_link(wc, ccre_df)
    output:
        out = temporary(config['ccre']['bed_gz'])

use rule gunzip as gunzip_ccre with:
    input:
        gz = config['ccre']['bed_gz']
    output:
        out = temporary(config['ccre']['bed'])

rule ccre_format:
    input:
        bed = config['ccre']['bed']
    resources:
        mem_gb = 10,
        threads = 1
    output:
        bed = config['ccre']['bed_format']
    run:
        df = pd.read_csv(input.bed, sep='\t',
                         header=None,
                         names=['Chromosome', 'Start','End'], usecols=[0,1,2])
        df['assay'] = wildcards.ccre_type
        df.to_csv(output.bed, index=False, sep='\t')

rule all_ccre:
    input:
        expand(rules.ccre_format.output,
               ccre_type=ccre_types,
               species=species)
