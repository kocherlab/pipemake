rule all:
    input:
        expand("MSF/AA/{sample}.fa", sample=config["samples"]),


rule translate_msf:
    input:
        "MSF/CDS/{sample}.fa",
    output:
        "MSF/AA/{sample}.fa",
    log:
        "logs/translate-msf/{sample}.log",
    singularity:
        "docker://aewebb/pipemake_utils:v1.4.3"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "translate-seq-file --input {input} --output {output} &> {log}"
