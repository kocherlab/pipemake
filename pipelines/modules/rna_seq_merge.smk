rule all:
    input:
        f"RNAseq/BAM/{config['species']}_{config['assembly_version']}.bam",


rule merge_bam:
    input:
        expand("RNAseq/BAM/Sorted/{sample}.sortedByCoord.bam", sample=config["samples"]),
    output:
        f"RNAseq/BAM/{config['species']}_{config['assembly_version']}.bam",
    log:
        f"logs/samtools/{config['species']}_{config['assembly_version']}.merge.log",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "samtools merge -@ {threads} -r {output} {input} 2> {log}"
