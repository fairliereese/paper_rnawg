species = ['human']
end_mode = ['tes']
pas_df = pd.read_csv(expand(config['pas']['encode_meta'],
                      species=species)[0],
                      sep='\t').set_index('File accession')

wildcard_constraints:
  encid='|'.join([re.escape(x) for x in pas_df.index.tolist()])

use rule dl_encid_2 as dl_pas with:
    params:
        file_ext = 'bam'
    output:
        out = config['pas']['bam']

rule pas_lapa_config:
    input:
        expand(config['pas']['bam'],
               encid=pas_df.index.tolist(),
               species=species)
    params:
        df_pas = pas_df
    threads: 1
    resources:
        mem_gb = 4
    output:
        config = config['pas']['config']
    script:
        "./lapa_pas_config.py"

rule lapa_call_ends_pas:
    input:
        config = config['pas']['config'],
        fa = config['ref']['talon']['fa'],
        gtf = config['ref']['lapa']['gtf'],
        chrom_sizes = config['ref']['talon']['chrom_sizes']
    params:
        replication_num_sample = 3,
        opref = config['pas']['ends'].rsplit('/', maxsplit=1)[0]+'/',
        lapa_cmd = lambda wc:get_lapa_settings(wc,
                                config['pas']['ends'],
                                'lapa_cmd'),
        lapa_end_temp = lambda wc:get_lapa_settings(wc,
                            config['pas']['ends'],
                            'temp_file')
    output:
        ends = config['pas']['ends']
    shell:
        """
        rm -rf {params.opref}
        {params.lapa_cmd} \
            --alignment {input.config} \
            --fasta {input.fa} \
            --annotation {input.gtf} \
            --chrom_sizes {input.chrom_sizes} \
            --replication_num_sample {params.replication_num_sample} \
            --output_dir {params.opref}
        if [ {params.lapa_end_temp} != {output.ends} ]
        then
            cp {params.lapa_end_temp} {output.ends}
        fi
        """

rule all_pas:
    input:
        expand(rules.lapa_call_ends_pas.output,
               species=species,
               end_mode=end_mode)
