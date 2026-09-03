rule all:
    input:
        expand("IsoSeq/BAM/Sorted/{sample}.sortedByCoord.bam", sample=config["samples"]),


rule sort_bam_isoseq:
    input:
        "IsoSeq/BAM/Aligned/{sample}.Aligned.bam",
    output:
        bam="IsoSeq/BAM/Sorted/{sample}.sortedByCoord.bam",
        index="IsoSeq/BAM/Sorted/{sample}.sortedByCoord.bam.bai",
    log:
        "logs/samtools/sort_{sample}.log",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=8000,
    threads: 4
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input} 2> {log}
        samtools index {output.bam} 2>> {log}
        """
