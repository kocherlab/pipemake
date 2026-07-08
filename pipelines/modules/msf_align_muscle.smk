rule all:
    input:
        expand("MSA/MUSCLE/{sample}.fa", sample=config["samples"]),


rule msf_align_muscle:
    input:
        "MSF/{sample}.fa",
    output:
        "MSA/MUSCLE/{sample}.fa",
    log:
        "logs/MUSCLE/{sample}.msf_align_muscle.log",
    singularity:
        "docker://aewebb/muscle:v5.3"
    resources:
        mem_mb=32000,
    threads: 1
    shell:
        "muscle -align {input} -output {output} 2> {log}"
