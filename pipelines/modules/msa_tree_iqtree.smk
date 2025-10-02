rule all:
    input:
        expand("iqtree/{sample}.fa.contree", sample=config["samples"]),


rule run_iqtree:
    input:
        "iqtree/{sample}.fa",
    output:
        "iqtree/{sample}.fa.contree",
    singularity:
        "docker://aewebb/iqtree:v3.0.1"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "iqtree -s {input} -alrt 1000 -B 1000 -T {threads}"
