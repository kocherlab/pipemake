rule all:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa.masked",


rule mask_assembly_counts_windowmasker:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        f"Assembly/WindowMasker/{config['species']}_{config['assembly_version']}.windowmasker.counts",
    resources:
        mem_mb=16000,
    threads: 1
    singularity:
        "docker://"
    shell:
        "windowmasker -mk_counts -in {input} -infmt fasta -sformat obinary -out {output}"


rule run_windowmasker:
    input:
        assembly=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
        counts=f"Assembly/WindowMasker/{config['species']}_{config['assembly_version']}.windowmasker.counts",
    output:
        f"Assembly/WindowMasker/{config['species']}_{config['assembly_version']}.fa.masked",
    resources:
        mem_mb=16000,
    threads: 1
    singularity:
        "docker://"
    shell:
        "windowmasker -in {input.assembly} -infmt fasta -ustat {input.counts} -dust T -outfmt fasta -out {output}"


rule store_windowmasker_assembly:
    input:
        f"Assembly/WindowMasker/{config['species']}_{config['assembly_version']}.fa.masked",
    output:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa.masked",
    localrule: True
    shell:
        "cp {input} {output}"
