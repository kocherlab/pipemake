rule all:
    input:
        expand("MSF/AA/{sample}.fasta", sample=config["samples"]),


rule translate_msf:
    input:
        "MSF/CDS/{sample}.fasta",
    output:
        "MSF/AA/{sample}.fasta",
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "translate-seq-file --input {input} --output {output}"
