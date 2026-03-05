rule all:
    input:
        expand("MSA/Trimmed/{sample}.fa", sample=config["samples"]),


rule trim_msa_clipkit:
    input:
        "MSA/Untrimmed/{sample}.fa",
    output:
        "MSA/Trimmed/{sample}.fa",
    params:
        method=config["clipkit_params"]["trim_method"],
        codon="--codon" if config["codon-msa"] else "",
    singularity:
        "docker://aewebb/clipkit:v2.11.4"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "clipkit {params.codon} --method {params.method} --output {output} {input}"
