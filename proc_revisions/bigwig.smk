def get_strand_flag(wc):
    """
    Get strandedness argument for bigwig creation
    """
    if wc.strand == 'fwd':
        flag = '--filterRNAstrand reverse'
    elif wc.strand == 'rev':
        flag = '--filterRNAstrand forward'
    return flag

rule bam_to_bw:
    resources:
        threads = 4,
        mem_gb = 64
    params:
        strand_flag = lambda wc:get_strand_flag(wc)
    conda:
        'deeptools'
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw} {params.strand_flag}
        """

use rule bam_to_bw as bw with:
    input:
        bam = config['lr']['talon']['bam_sort'],
        bai = config['lr']['talon']['bam_sort']
    output:
        bw = config['lr']['bw']

rule all_bigwig:
    input:
        expand(expand(config['lr']['bw'],
                      zip,
                      species=lr_df.species.tolist(),
                      dataset=lr_df.dataset.tolist(),
                      allow_missing=True),
                      strand=['fwd', 'rev'])
