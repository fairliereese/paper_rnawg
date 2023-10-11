rule talon_init:
	resources:
		mem_gb = 32,
		threads = 16
	shell:
		'talon_initialize_database \
    		--f {input.gtf} \
    		--g {params.genome_ver} \
    		--a {params.annot_ver} \
    		--l 0 \
    		--idprefix TALON \
    		--5p 500 \
    		--3p 300 \
    		--o {params.opref}'

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
            -v 2 > {output.debug_log}
        mv ${{ref_db}} {params.opref}_talon.db
        """
