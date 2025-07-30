rule all:
    input:
        expand("reSEQ/BAM/Sorted/{sample}.sortedByCoord.bam", sample=config["samples"]),
        expand(
            "reSEQ/BAM/Sorted/{sample}.sortedByCoord.bam.bai", sample=config["samples"]
        ),


rule sort_bam_reseq:
    input:
        "reSEQ/BAM/Aligned/{sample}.Aligned.bam",
    output:
        bam="reSEQ/BAM/Sorted/{sample}.sortedByCoord.bam",
        index="reSEQ/BAM/Sorted/{sample}.sortedByCoord.bam.bai",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=8000,
    threads: 4
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input}
        samtools index {output.bam}
        """
