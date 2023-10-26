species = ['human', 'mouse']

use rule dl as dl_fantom_cage with:
    params:
        link = lambda wc:config['fantom'][wc.species]['link']
    output:
        out = temporary(config['fantom']['bed_old_gz'])

use rule gunzip as gz_fantom_cage with:
    input:
        gz = config['fantom']['bed_old_gz']
    output:
        out = temporary(config['fantom']['bed_old'])

rule liftover_fantom:
    input:
        bed = config['fantom']['bed_old']
    resources:
        threads = 1,
        mem_gb = 16
    log: config['fantom']['log']
    output:
        bed = config['fantom']['bed']
    script:
        "./liftover_fantom.py"

rule all_fantom:
    input:
        expand(rules.liftover_fantom.output,
               species=species)
