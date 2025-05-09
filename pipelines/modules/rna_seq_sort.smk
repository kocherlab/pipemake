rule all:
    input:
        expand("RNAseq/BAM/Sorted/{sample}.sortedByCoord.bam", sample=config["samples"]),


rule sort_bam_rnaseq:
    input:
        "RNAseq/BAM/Aligned/{sample}.Aligned.bam",
    output:
        "RNAseq/BAM/Sorted/{sample}.sortedByCoord.bam",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=8000,
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output} {input}"
