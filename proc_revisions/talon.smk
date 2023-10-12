from sm_utils import *

############# config stuff
configfile: 'config.yml'

datasets_per_talon_run = config['params']['talon']['datasets_per_run']
min_transcript_len = config['params']['talon']['min_transcript_len']
int_priming_a_range = config['params']['talon']['int_priming_a_range']
max_5_dist = config['params']['talon']['max_5_dist']
max_3_dist = config['params']['talon']['max_3_dist']
max_frac_a = config['params']['talon']['max_frac_a']
min_count = config['params']['talon']['min_count']
min_datasets = config['params']['talon']['min_datasets']


lr_df = process_lr_metadata(config['lr']['meta'],
                            ['human', 'mouse'],
                            datasets_per_talon_run)

############# snakemake settings stuff
ruleorder:
    talon_first > talon_seq

def get_talon_run_file(talon_run, cfg_entry):
    """
    Get the config entry associated with the talon
    run number. Goal of this function is
    to make my expand statements more reproducible
    """
    files = expand(cfg_entry,
                   zip,
                   talon_run=int(talon_run),
                   allow_missing=True)
    assert len(files) == 1
    return files[0]

def get_talon_run_info(wc, df, cfg_entry, dataframe=False):
    """
    Get all files for a talon run

    Parameters:
        dataframe (bool): False if it should return just the
            list of input files for this TALON run,
            True if it should return the whole DF
    """
    temp = df.copy(deep=True)
    temp = temp.loc[(temp.species==wc.species)&\
                  (temp.talon_run==int(wc.talon_run))]
    datasets = temp.dataset.tolist()
    species = temp.species.tolist()
    files = expand(cfg_entry,
                   zip,
                   dataset=datasets,
                   species=species)
    temp['talon_file'] = files
    if not dataframe:
        return files
    else:
        return temp

################################################################################
################################ Template rules ################################
################################################################################

rule talon_init:
	resources:
		mem_gb = 32,
		threads = 16
	shell:
		"""
		talon_initialize_database \
    		--f {input.gtf} \
    		--g {params.genome_ver} \
    		--a {params.annot_ver} \
    		--l {params.min_transcript_len} \
    		--idprefix TALON \
    		--5p {params.max_5_dist} \
    		--3p {params.max_3_dist} \
    		--o {params.opref}
		"""

rule talon:
    resources:
        mem_gb = 256,
        threads = 30
    shell:
        """
        ref_db={input.ref}_{wildcards.species}
        cp {input.ref} ${{ref_db}}
        talon \
            --f {input.config} \
            --db ${{ref_db}} \
            --build {params.genome_ver} \
            --tmpDir {params.opref}_temp/ \
            --threads {resources.threads} \
            --create_novel_spliced_genes \
            --o {params.opref} \
            -v 1 > {output.debug_log}
        mv ${{ref_db}} {params.opref}_talon.db
        """

rule talon_annot:
    resources:
        mem_gb = 128,
        threads = 1
    shell:
        """
		talon_fetch_reads \
            --db {input.db} \
            --build {params.genome_ver} \
            --o {params.opref}
		"""

rule talon_abundance:
    resources:
        mem_gb = 128,
        threads = 1
    shell:
        """
		talon_abundance \
            --db {input.db} \
            -a {params.annot_ver} \
            -b {params.genome_ver} \
            --o {params.opref}
		"""

rule talon_filter:
    resources:
        mem_gb = 128,
        threads = 1
    shell:
        """
		talon_filter_transcripts \
            --db {input.db} \
            -a {params.annot_ver} \
            --maxFracA {params.max_frac_a}\
            --minCount {params.min_count} \
            --minDatasets {params.min_datasets} \
            --o {output.pass_list}
		"""

rule talon_filtered_abundance:
    resources:
        mem_gb = 128,
        threads = 1
    shell:
        """
		talon_abundance \
            --db {input.db} \
            -a {params.annot_ver} \
            -b {params.genome_ver} \
            --whitelist {input.pass_list} \
            --o {params.opref}
		"""

rule talon_gtf:
    resources:
        mem_gb = 128,
        threads = 1
    shell:
        """
		talon_create_GTF \
            --db {input.db} \
            -a {params.annot_ver} \
            -b {params.genome_ver} \
            --whitelist {input.pass_list} \
            --o {params.opref}
		"""

################################################################################
#################################### TALON #####################################
################################################################################

# label reads and create indexed BAM files to retain
rule talon_label:
    input:
        bam = config['lr']['bam'],
        fa = config['ref']['talon']['fa']
    resources:
        mem_gb = 64,
        threads = 1
    params:
        opref = config['lr']['talon']['sam_label'].rsplit('_labeled.sam')[0],
        int_priming_a_range = int_priming_a_range
    output:
        bam = temporary(config['lr']['talon']['sam_label']),
        tsv = temporary(config['lr']['talon']['read_labels'])
    shell:
        """
        talon_label_reads \
            --f {input.bam} \
            --g {input.fa} \
            --tmpDir {params.opref}_tmp/ \
            --ar {params.int_priming_a_range}  \
            --deleteTmp  \
            --o {params.opref}
        """

use rule sam_to_bam as lr_sam_to_bam with:
    input:
        sam = config['lr']['talon']['sam_label']
    output:
        bam = temporary(config['lr']['talon']['bam'])

use rule sort_bam as lr_sort_bam with:
    input:
        bam = config['lr']['talon']['bam']
    output:
        bam = config['lr']['talon']['bam_sort']

use rule index_bam as lr_index_bam with:
    input:
        bam = config['lr']['talon']['bam_sort']
    output:
        ind = config['lr']['talon']['bam_ind']

# initialize TALON db
use rule talon_init as talon_init_db with:
    input:
        gtf = config['ref']['talon']['gtf']
    params:
        min_transcript_len = min_transcript_len,
        max_5_dist = max_5_dist,
        max_3_dist = max_3_dist,
        genome_ver = lambda wc:config['ref'][wc.species]['fa_ver'],
        annot_ver = lambda wc:config['ref'][wc.species]['gtf_ver'],
        opref = config['lr']['talon']['ref_db'].rsplit('.db', maxsplit=1)[0]
    output:
        db = config['lr']['talon']['ref_db']

# get config files for each talon run
rule talon_config:
    input:
        files = lambda wc:get_talon_run_info(wc, lr_df,
                            cfg_entry=config['lr']['talon']['bam_sort'],
                            dataframe=False)
    resources:
        threads = 1,
        mem_gb = 1
    params:
        df = lr_df
    output:
        config = config['lr']['talon']['config']
    run:
        cfg_df = get_talon_run_info(wildcards, params.df,
                                    config['lr']['talon']['bam_sort'],
                                    dataframe=True)
        cfg_df = cfg_df[['dataset', 'sample', 'platform', 'talon_file']]
        cfg_df.to_csv(output.config, header=None, sep=',', index=False)


# first talon run
use rule talon as talon_first with:
    input:
        ref = config['lr']['talon']['ref_db'],
        config = get_talon_run_file(0, config['lr']['talon']['config'])
    params:
        genome_ver = lambda wc:config['ref'][wc.species]['fa_ver'],
        opref = get_talon_run_file(0, config['lr']['talon']['db']).rsplit('_talon.db',
                                   maxsplit=1)[0]
    output:
        db = get_talon_run_file(0, config['lr']['talon']['db']),
        read_annot = get_talon_run_file(0, config['lr']['talon']['annot']),
        debug_log = get_talon_run_file(0, config['lr']['talon']['debug_log'])

# sequential talon runs
use rule talon as talon_seq with:
    input:
        ref = lambda wc:get_talon_run_file(int(wc.talon_run)-1,
                        config['lr']['talon']['db']),
        config = config['lr']['talon']['config']
    params:
        genome_ver = lambda wc:config['ref'][wc.species]['fa_ver'],
        opref = config['lr']['talon']['db'].rsplit('_talon.db', maxsplit=1)[0]
    output:
        db = config['lr']['talon']['db'],
        read_annot = config['lr']['talon']['annot'],
        debug_log = config['lr']['talon']['debug_log']

################################################################################
################################# TALON output #################################
################################################################################

def get_max_talon_run_cfg(wc, df, cfg_entry):
	"""
	Get the last TALON run given the wildcards and
	output the corresponding file needed
	"""
	temp = df.loc[df.species==wc.species].copy(deep=True)
	max_talon_run = temp.max_talon_run.tolist()[0]
	file = expand(cfg_entry,
				  zip,
				  talon_run=max_talon_run,
				  species=wc.species)
	assert len(file) == 1
	return file

use rule talon_annot as talon_annot_full with:
    input:
        db = lambda wc:get_max_talon_run_cfg(wc, lr_df,config['lr']['talon']['db'])
    params:
        opref = config['lr']['talon']['full_annot'].rsplit('_talon', maxsplit=1)[0],
        genome_ver = lambda wc:config['ref'][wc.species]['fa_ver']
    output:
        read_annot = config['lr']['talon']['full_annot']

use rule talon_abundance as talon_ab_full with:
    input:
        db = lambda wc:get_max_talon_run_cfg(wc, lr_df, config['lr']['talon']['db'])
    params:
        genome_ver = lambda wc:config['ref'][wc.species]['fa_ver'],
        annot_ver = lambda wc:config['ref'][wc.species]['gtf_ver'],
        opref = config['lr']['talon']['ab'].rsplit('_talon', maxsplit=1)[0]
    output:
        ab = config['lr']['talon']['ab']

use rule talon_filter as talon_filt_full with:
    input:
        db = lambda wc:get_max_talon_run_cfg(wc, lr_df, config['lr']['talon']['db'])
    params:
        annot_ver = lambda wc:config['ref'][wc.species]['gtf_ver'],
        max_frac_a = max_frac_a,
        min_count = min_count,
        min_datasets = min_datasets
    output:
        pass_list = config['lr']['talon']['pass_list']

use rule talon_filtered_abundance as talon_filt_ab_full with:
    input:
        db = lambda wc:get_max_talon_run_cfg(wc, lr_df, config['lr']['talon']['db']),
        pass_list = config['lr']['talon']['pass_list']
    params:
        opref = config['lr']['talon']['filt_ab'].rsplit('_talon', maxsplit=1)[0],
        genome_ver = lambda wc:config['ref'][wc.species]['fa_ver'],
        annot_ver = lambda wc:config['ref'][wc.species]['gtf_ver']
    output:
        ab = config['lr']['talon']['filt_ab']

use rule talon_gtf as talon_gtf_full with:
    input:
        db = lambda wc:get_max_talon_run_cfg(wc, lr_df, config['lr']['talon']['db']),
        pass_list = config['lr']['talon']['pass_list']
    params:
        opref = config['lr']['talon']['gtf'].rsplit('_talon', maxsplit=1)[0],
        genome_ver = lambda wc:config['ref'][wc.species]['fa_ver'],
        annot_ver = lambda wc:config['ref'][wc.species]['gtf_ver']
    output:
        gtf = config['lr']['talon']['gtf']
