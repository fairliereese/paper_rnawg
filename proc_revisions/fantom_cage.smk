species = ['human']

use rule dl as dl_fantom_cage with:
    params:
        link = config['fantom']['link']
    output:
        out = temporary(config['fantom']['bed_hg19_gz'])

use rule gunzip as gz_fantom_cage with:
    input:
        gz = config['fantom']['bed_hg19_gz']
    output:
        out = config['fantom']['bed_hg19']

rule liftover_fantom:
    input:
        bed = config['fantom']['bed_hg19']
    resources:
        threads = 1,
        mem_gb = 16
    output:
        bed = config['fantom']['bed']
    script:
        "./liftover_fantom.py"

rule all_fantom:
    input:
        expand(rules.liftover_fantom.output,
               species=species)
