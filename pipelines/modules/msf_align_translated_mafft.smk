rule all:
    input:
        expand("MSA/AA/MAFFT/{sample}.fasta", sample=config["samples"]),


rule align_msa_mafft:
    input:
        unaligned_aa="MSF/AA/{sample}.fasta",
    output:
        aligned_aa="MSA/AA/MAFFT/{sample}.fasta",
    singularity:
        "docker://aewebb/mafft:v7.525"
    resources:
        mem_mb=8000,
    threads: 4
    shell:
        "mafft-linsi --thread {threads} {input.unaligned_aa} > {output.aligned_aa}"
