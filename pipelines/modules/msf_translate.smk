rule all:
    input:
        expand("MSF/AA/{sample}.fa", sample=config["samples"]),


rule translate_msf:
    input:
        "MSF/CDS/{sample}.fa",
    output:
        "MSF/AA/{sample}.fa",
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.8"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "translate-seq-file --input {input} --output {output}"
