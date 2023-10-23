species = ['human']

use rule dl as dl_polyasite_atlas with:
    params:
        link = config['polya_atlas']['link']
    output:
        out = temporary(config['polya_atlas']['bed_gz'])

use rule gunzip as gz_polyasite_atlas with:
    input:
        gz = config['polya_atlas']['bed_gz']
    output:
        out = config['polya_atlas']['bed']

rule all_polyasite:
    input:
        expand(config['polya_atlas']['bed'],
               species=species)
