rule all:
    input:
        f"IsoSeq/BAM/{config['species']}_{config['assembly_version']}.bam",


rule merge_bam:
    input:
        expand("IsoSeq/BAM/Sorted/{sample}.sortedByCoord.bam", sample=config["samples"]),
    output:
        f"IsoSeq/BAM/{config['species']}_{config['assembly_version']}.bam",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "samtools merge -@ {threads} -r {output} {input}"
