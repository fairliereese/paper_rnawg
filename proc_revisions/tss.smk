# cage stuff
species = ['human']
cage_df = pd.read_csv(expand(config['cage']['encode_meta'],
                      species=species)[0],
                      sep='\t').set_index('File accession')

rampage_df = pd.read_csv(expand(config['rampage']['encode_meta'],
                    species=species)[0],
                    sep='\t').set_index('File accession')

use rule dl_encid_gz_2 as dl_cage with:
    params:
        file_ext = 'bed'
    output:
        out = config['cage']['bed_gz']

rule merge_cage:
    input:
        beds = expand(config['cage']['bed_gz'],
                      encid=cage_df.index.tolist(),
                      species=species)
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
        out = config['rampage']['bed_gz']

rule merge_rampage:
    input:
        beds = expand(config['rampage']['bed_gz'],
                      encid=rampage_df.index.tolist(),
                      species=species)
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
               species=species),
        expand(rules.merge_rampage.output,
               species=species)
