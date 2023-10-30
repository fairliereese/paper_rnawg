import swan_vis as swan
from utils import *

################################################################################
################################## MANE ########################################
################################################################################
rule get_exp_genes:
    input:
        ab = rules.swan_init.input.ab
    resources:
        mem_gb = 16,
        threads = 8
    params:
        min_tpm = 1,
        gene_subset = 'polya',
    output:
        ofile = config['lr']['mane']['exp_gene_subset']
    run:
        get_exp_gene_subset(input.ab,
                            params.min_tpm,
                            wildcards.obs_col,
                            output.ofile,
                            wildcards.species)

rule get_pi_tpm:
    input:
        swan_file = rules.swan_init.output.sg,
        enc_meta = config['lr']['encode_meta']
    resources:
        mem_gb = 32,
        threads = 8
    params:
        odir = config['lr']['mane']['pi_tpm']['tss'].rsplit('/', maxsplit=1)[0]
    output:
        config['lr']['mane']['pi_tpm']['tss'],
        config['lr']['mane']['pi_tpm']['ic'],
        config['lr']['mane']['pi_tpm']['tes'],
        config['lr']['mane']['pi_tpm']['triplet']
    run:
        sg = swan.read(input.swan_file)

        # add biosample info in case we're using that
        meta_df = pd.read_csv(input.enc_meta, sep='\t')
        meta_df = format_metadata_col(meta_df, 'Biosample term name', 'biosample')
        meta_df = meta_df[['Experiment accession', 'biosample']]
        meta_df.drop_duplicates(inplace=True)
        meta_df.rename({'Experiment accession': 'ENCODE_experiment_id'},
                       axis=1, inplace=True)
        sg.adata.obs = sg.adata.obs.merge(meta_df,
                  how='left',
                  on='ENCODE_experiment_id')
        sg.tss_adata.obs = sg.tss_adata.obs.merge(meta_df,
                how='left',
                on='ENCODE_experiment_id')
        sg.ic_adata.obs = sg.ic_adata.obs.merge(meta_df,
                how='left',
                on='ENCODE_experiment_id')
        sg.tes_adata.obs = sg.tes_adata.obs.merge(meta_df,
                how='left',
                on='ENCODE_experiment_id')

        # get the actual table
        get_pi_tpm_tables(sg, wildcards.obs_col, params.odir)


rule all_mane:
    input:
        config['lr']['mane']['exp_gene_subset'],
        config['lr']['mane']['pi_tpm']['tss']
