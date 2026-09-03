rule all:
    input:
        expand("Alignment/miniprot/AA/{sample}.fa", sample=config["samples"]),


rule fasta_translate_for_miniprot:
    input:
        "Alignment/Query/CDS/{sample}.fa",
    output:
        "Alignment/Query/AA/{sample}.fa",
    log:
        "logs/translate-seq-file/{sample}.log",
    singularity:
        "docker://aewebb/pipemake_utils:v1.4.3"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "translate-seq-file --input {input} --output {output} &> {log}"


rule index_assembly_for_miniprot:
    input:
        assembly_fasta=config["assembly_input"],
    output:
        temp("Index/miniprot/Assembly.mmi"),
    log:
        "logs/miniprot/index-assembly.log",
    singularity:
        "docker://aewebb/miniprot:v0.18-r281"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "miniprot -t {threads} -d {output} {input} &> {log}"


rule fasta_align_miniprot:
    input:
        protein_fasta="Alignment/Query/AA/{sample}.fa",
        assembly_fasta="Index/miniprot/Assembly.mmi",
    output:
        "Alignment/miniprot/GFF/{sample}.gff",
    log:
        "logs/miniprot/{sample}.align.log",
    params:
        species=config["species"],
    singularity:
        "docker://aewebb/miniprot:v0.18-r281"
    resources:
        mem_mb=8000,
    threads: 2
    shell:
        "miniprot -t {threads} -P {params.species}-{wildcards.sample}- --gff-only {input.assembly_fasta} {input.protein_fasta} > {output} 2> {log}"


rule miniprot_gff_to_proteins:
    input:
        gff_file="Alignment/miniprot/GFF/{sample}.gff",
        assembly_fasta=config["assembly_input"],
        assembly_idx=f"{config['assembly_input']}.fai",
    output:
        "Alignment/miniprot/AA/{sample}.fa",
    log:
        "logs/gffread/{sample}.proteins.log",
    singularity:
        "docker://aewebb/gffread:v0.12.7"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gffread -V -y {output} -g {input.assembly_fasta} {input.gff_file} &> {log}"
