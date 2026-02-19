rule all:
    input:
        expand("MSA/AA/{sample}.fa", sample=config["samples"]),


rule msf_align_muscle:
    input:
        "MSF/AA/{sample}.fa",
    output:
        "MSA/AA/{sample}.fa",
    log:
        "logs/MUSCLE/{sample}.msf_align_muscle.log",
    singularity:
        "docker://aewebb/muscle:v5.3"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "muscle -align {input} -output {output}"
