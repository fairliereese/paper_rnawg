rule sort_bam:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools sort \
            --threads {resources.threads} \
            -O bam {input.bam} > {output.bam}
        """

rule index_bam:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools index -@ {resources.threads} {input.bam}
        """
