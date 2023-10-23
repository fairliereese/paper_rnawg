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

rule dl_encid_2:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "wget https://www.encodeproject.org/files/{wildcards.encid}/@@download/{wildcards.encid}.{params.file_ext} -O {output.out}"


# download w/o an alias for the name
rule dl_encid_gz_2:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "wget https://www.encodeproject.org/files/{wildcards.encid}/@@download/{wildcards.encid}.{params.file_ext}.gz -O {output.out}"
