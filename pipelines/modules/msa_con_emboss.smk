rule all:
    input:
        expand("Consensus/{sample}.fa", sample=config["samples"]),


rule concencus_sequence_emboss:
    input:
        "MSA/{sample}.fa",
    output:
        "Consensus/{sample}.fa",
    params:
        plurality=config["emboss_cons_params"]["plurality"],
    singularity:
        "docker://biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "/usr/lib/emboss/cons -plurality {params.plurality} -sequence {input} -outseq {output} -name {wildcards.sample}"
