import itertools

wildcard_constraints:
    obs_cond1='|'.join([re.escape(x) for x in lr_df['sample'].unique().tolist()]),
    obs_cond2='|'.join([re.escape(x) for x in lr_df['sample'].unique().tolist()]),
    feat='|'.join([re.escape(x) for x in ['tss', 'tes', 'ic', 'iso']]),
    obs_col='sample|dataset'

def get_tc_mouse_samples(config):
    m_lib_meta = expand(config['lr']['meta'], species='mouse')[0]
    temp = pd.read_csv(m_lib_meta, sep='\t')
    tc_tissues = ['muscle', 'hippocampus', 'cortex', 'adrenal gland', 'heart']
    tc_times = ['18-20mo', '2mo', '4d', '25d', '14d', '36d', '10d']
    s = temp.loc[(temp.general_tissue_cell_type.isin(tc_tissues))&\
             (temp.age.isin(tc_times)), 'sample'].unique().tolist()
    return s

def get_du_tc_cfg_entries():
    """
    Get the cfg entries for running du tests
    """
    obs_col = 'sample'
    feats = ['tss', 'tes', 'ic', 'iso']
    s = get_tc_mouse_samples(config)
    combos = [c for c in itertools.combinations(s, 2) if c[0].split('_')[0]==c[1].split('_')[0]]
    obs_cond1 = [c[0] for c in combos]
    obs_cond2 = [c[1] for c in combos]
    files = expand(expand(config['lr']['analysis']['du'],
                  zip,
                  obs_cond1=obs_cond1,
                  obs_cond2=obs_cond2,
                  allow_missing=True),
                  obs_col='sample',
                  species='mouse',
                  feat=feats)
    return files

obs_col = 'sample'
feats = ['tss', 'tes', 'ic', 'iso']
s = get_tc_mouse_samples(config)
combos = [c for c in itertools.combinations(s, 2) if c[0].split('_')[0]==c[1].split('_')[0]]
obs_cond1 = [c[0] for c in combos]
obs_cond2 = [c[1] for c in combos]
files = expand(expand(config['lr']['analysis']['du'],
              zip,
              obs_cond1=obs_cond1,
              obs_cond2=obs_cond2,
              allow_missing=True),
              obs_col='sample',
              species='mouse',
              feat=feats)

rule swan_die:
    input:
        sg = expand(config['lr']['swan']['sg'], species='mouse')[0],
        meta = expand(config['lr']['meta'], species='mouse')[0]
        # sg = config['lr']['swan']['sg'],
        # meta = config['lr']['meta'],
    resources:
        mem_gb = 128,
        threads = 8
    output:
        out = expand(config['lr']['analysis']['du'],
                     species='mouse')
    run:
        sg = swan.read(input.sg)
        die, genes = sg.die_gene_test(obs_col=wildcards.obs_col,
                                      obs_conditions=[wildcards.obs_cond1,
                                                      wildcards.obs_cond2],
                                      kind=wildcards.feat)
        die.to_csv(output.out, sep='\t')



rule all_du:
    input:
        expand(expand(config['lr']['analysis']['du'],
                      zip,
                      obs_cond1='adrenal_10d',
                      obs_cond2='adrenal_14d',
                      allow_missing=True),
                      obs_col='sample',
                      species='mouse',
                      feat='iso')
        # expand(config['lr']['analysis']['du'],
        #               zip,
        #               obs_cond1='adrenal_10d',
        #               obs_cond2='adrenal_14d',
        #               obs_col='sample',
        #               species='mouse',
        #               feat='iso')
