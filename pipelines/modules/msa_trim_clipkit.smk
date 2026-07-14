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
        codon="--codon" if "codon_msa" in config else "",
    singularity:
        "docker://aewebb/clipkit:v2.11.4"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "clipkit {params.codon} --mode {params.method} --output {output} {input}"
