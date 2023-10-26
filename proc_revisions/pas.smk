species = ['human', 'mouse']
end_mode = ['tes']

pas_df = pd.read_csv(expand(config['pas']['encode_meta'],
                      species='human')[0],
                      sep='\t').set_index('File accession')
pas_df['species'] = 'human'
temp = pd.read_csv(expand(config['pas']['encode_meta'],
                      species='mouse')[0],
                      sep='\t').set_index('File accession')
temp['species'] = 'mouse'
pas_df = pd.concat([pas_df, temp], axis=0)

def get_species_files(wc, df, col=None):
    temp = df.loc[df.species==wc.species].copy(deep=True)
    if col:
        vals = temp[col].tolist()
        return vals
    else:
        return temp

wildcard_constraints:
  encid='|'.join([re.escape(x) for x in pas_df.index.tolist()])

use rule dl_encid_2 as dl_pas with:
    params:
        file_ext = 'bam'
    output:
        out = config['pas']['bam']

rule pas_lapa_config:
    input:
        lambda wc:expand(config['pas']['bam'],
               encid=get_species_files(wc, pas_df).index.tolist(),
               species=wc.species)
    params:
        df_pas = pas_df
    resources:
        threads = 1,
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

rule format_lapa_ends:
    input:
        bed = config['pas']['ends']
    resources:
        mem_gb = 4,
        threads = 1
    output:
        bed_formatted = config['pas']['ends_formatted']
    run:
        format_lapa_ends(input.bed, output.bed_formatted)


rule all_pas:
    input:
        expand(rules.format_lapa_ends.output,
               species='human',
               end_mode='tss')
        # expand(rules.format_lapa_ends.output,
        #        species=species,
        #        end_mode=end_mode)
