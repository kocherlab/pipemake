rule all:
    input:
        expand("HiFi/PAF/{sample}.paf", sample=config["samples"]),


rule hifi_assembly_align_minimap2:
    input:
        hifi_assembly_fasta="HiFi/Assembly/{sample}.fa",
        assembly_fasta=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        "HiFi/PAF/{sample}.paf",
    singularity:
        "docker://aewebb/minimap2:v2.28"
    resources:
        mem_mb=32000,
    threads: 16
    shell:
        "minimap2 -x asm5 -t {threads} {input.assembly_fasta} {input.hifi_assembly_fasta} > {output}"
