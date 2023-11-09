lrgasp_species = ['human']

meta = pd.read_csv(expand(config['lrgasp_cage']['metadata'],
                   species=lrgasp_species)[0],
                   sep='\t')

wildcard_constraints:
   lrgasp_cage_dataset='|'.join([re.escape(x) for x in meta.dataset.tolist()])


def get_col_from_dataset(dataset, meta, col):
    return meta.loc[meta.dataset==dataset, col].values[0]

use rule dl as dl_lrgasp_cage with:
    params:
        link = lambda wc:get_col_from_dataset(wc.lrgasp_cage_dataset,
                                             meta,
                                             'link')
    output:
        out = temporary(config['lrgasp_cage']['bed_gz'])

use rule gunzip as gz_lrgasp_cage with:
    input:
        gz = config['lrgasp_cage']['bed_gz']
    output:
        out = temporary(config['lrgasp_cage']['bed'])

rule merge_lrgasp_cage:
    input:
        beds = expand(config['lrgasp_cage']['bed'],
                      lrgasp_cage_dataset=meta.dataset.tolist(),
                      species=lrgasp_species)
    resources:
        mem_gb = 8,
        threads = 1
    output:
        bed = config['lrgasp_cage']['merged']
    script:
        "merge_beds.py"

rule all_lrgasp_cage:
    input:
        expand(config['lrgasp_cage']['merged'],
               species=lrgasp_species)
