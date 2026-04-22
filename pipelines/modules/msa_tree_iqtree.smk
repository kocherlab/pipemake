rule all:
    input:
        expand("iqtree/{sample}.fa.contree", sample=config["samples"]),


rule reconstruct_tree_iqtree:
    input:
        "iqtree/{sample}.fa",
    output:
        "iqtree/{sample}.fa.contree",
    params:
        msa=subpath(output[0], strip_suffix=".contree"),
    singularity:
        "docker://aewebb/iqtree:v3.0.1"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        """
        if [ "{input}" != "{params.msa}" ]; then
            cp {input} {params.msa}
        fi
        iqtree -s {params.msa} -alrt 1000 -B 1000 -T {threads}
        """
