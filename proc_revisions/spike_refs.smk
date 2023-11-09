
################################################################################
######################### Spike download / processing ##########################
################################################################################

use rule dl as dl_sirv_gtf with:
  params:
    link = config['ref']['spike']['sirv_gtf_link']
  output:
    out = temporary(config['ref']['spike']['sirv_gtf_gz'])

use rule gunzip as gz_sirv_gtf with:
  input:
    gz = config['ref']['spike']['sirv_gtf_gz']
  output:
    out = config['ref']['spike']['sirv_gtf']

use rule dl as dl_sirv_fa with:
  params:
    link = config['ref']['spike']['sirv_fa_link']
  output:
    out = temporary(config['ref']['spike']['sirv_fa_gz'])

use rule gunzip as gz_sirv_fa with:
  input:
    gz = config['ref']['spike']['sirv_fa_gz']
  output:
    out = config['ref']['spike']['sirv_fa']

use rule dl as dl_ercc_fa with:
  params:
    link = config['ref']['spike']['ercc_fa_link']
  output:
    out = temporary(config['ref']['spike']['ercc_fa_gz'])

use rule gunzip as gz_ercc_fa with:
  input:
    gz = config['ref']['spike']['ercc_fa_gz']
  output:
    out = config['ref']['spike']['ercc_fa']

rule mkref_ercc_gtf:
    input:
        fa = config['ref']['spike']['ercc_fa']
    resources:
        mem_gb = 16,
        threads = 1
    output:
        gtf = config['ref']['spike']['ercc_gtf']
    shell:
        """
        python ../scripts/merge_encode_annotations.py \
            -o temp_reformat_ercc.gtf \
            {input.fa}

        python ../scripts/reformat_gtf.py \
            -o {output.gtf} \
            -gtf temp_reformat_ercc.gtf

        rm temp_reformat_ercc.gtf
        """

# use rule get_chrom_sizes as get_spike_chrom_sizes with:
#     input:
#         fa = config['ref']['spike']['spike_fa']
#     output:
#         chrom_sizes = config['ref']['spike']['chrom_sizes']
#
# rule get_spike_chrom_bed:
#     input:
#         chrom_sizes = config['ref']['spike']['chrom_sizes']
#     resources:
#         mem_gb = 16,
#         threads = 1
#     output:
#         bed = config['ref']['spike']['chrom_bed']
#     run:
#         df = pd.read_csv(input.chrom_sizes,
#                          sep='\t',
#                          header=None,
#                          names=['chr', 'stop'])
#         df['start'] = 0
#         df = df[['chr', 'start', 'stop']]
#         df.to_csv(output.bed,
#                   sep='\t',
#                   header=None,
#                   index=False)
