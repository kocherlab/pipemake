rule all:
    input:
        expand("Alignment/miniprot/GFF/{sample}.gff", sample=config["samples"]),


rule fasta_translate_for_miniprot:
    input:
        "Alignment/Query/CDS/{sample}.fa",
    output:
        "Alignment/Query/AA/{sample}.fa",
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.8"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "translate-seq-file --input {input} --output {output}"

rule fasta_align_miniprot:
    input:
        protein_fasta="Alignment/Query/AA/{sample}.fa",
        assembly_fasta=config["assembly_input"],
    output:
        "Alignment/miniprot/GFF/{sample}.gff"
    log:
        "logs/miniprot/{sample}.log"
    singularity:
        "docker://aewebb/miniprot:v0.18-r281"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        'miniprot -t {threads} --gff-only {input.assembly_fasta} {input.protein_fasta} > {output} 2> {log}'
