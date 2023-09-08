################################################################################
########################### Rule definitions ###################################
################################################################################

rule dl:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "wget -O {output.out} {params.link}"

rule gunzip:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "gunzip -c {input.gz} > {output.out}"

rule dl_encid_gz:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.{params.file_ext}.gz -O {output.out}"

rule dl_encid:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "wget https://www.encodeproject.org/files/{params.encid}/@@download/{params.encid}.{params.file_ext} -O {output.out}"
