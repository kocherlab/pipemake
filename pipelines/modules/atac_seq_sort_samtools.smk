rule all:
    input:
        expand(
            "ATAC_seq/BAM/Sorted/{sample}.sortedByCoord.bam", sample=config["samples"]
        ),


rule sort_bam:
    input:
        "ATAC_seq/BAM/Aligned/{sample}.Aligned.bam",
    output:
        temp("ATAC_seq/BAM/Sorted/{sample}.sortedByCoord.bam"),
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output} {input}"
