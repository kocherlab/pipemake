rule all:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa.fai",


rule index_assembly:
    input:
        config["assembly_input"],
    output:
        f"{config['assembly_input']}.fai",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "samtools faidx {input}"
