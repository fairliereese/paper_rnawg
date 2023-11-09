# cage stuff
tss_species = ['human']
cage_df = pd.read_csv(expand(config['cage']['encode_meta'],
                      species=tss_species)[0],
                      sep='\t').set_index('File accession')

rampage_df = pd.read_csv(expand(config['rampage']['encode_meta'],
                    species=tss_species)[0],
                    sep='\t').set_index('File accession')

wildcard_constraints:
    encid='|'.join([re.escape(x) for x in cage_df.index.tolist()+rampage_df.index.tolist()]),

use rule dl_encid_gz_2 as dl_cage with:
    params:
        file_ext = 'bed'
    output:
        out = temporary(config['cage']['bed_gz'])

use rule gunzip as gz_cage with:
    input:
        gz = config['cage']['bed_gz']
    output:
        out = temporary(config['cage']['bed'])

rule merge_cage:
    input:
        beds = expand(config['cage']['bed'],
                      encid=cage_df.index.tolist(),
                      species=tss_species)
    resources:
        mem_gb = 8,
        threads = 1
    output:
        bed = config['cage']['merged']
    script:
        "merge_beds.py"

use rule dl_encid_gz_2 as dl_rampage with:
    params:
        file_ext = 'bed'
    output:
        out = temporary(config['rampage']['bed_gz'])

use rule gunzip as gz_rampage with:
    input:
        gz = config['rampage']['bed_gz']
    output:
        out = temporary(config['rampage']['bed'])

rule merge_rampage:
    input:
        beds = expand(config['rampage']['bed'],
                      encid=rampage_df.index.tolist(),
                      species=tss_species)
    resources:
        threads = 1,
        mem_gb = 8
    output:
        bed = config['rampage']['merged']
    script:
        "merge_beds.py"

rule all_tss:
    input:
        expand(rules.merge_cage.output,
               species=tss_species),
        expand(rules.merge_rampage.output,
               species=tss_species)
