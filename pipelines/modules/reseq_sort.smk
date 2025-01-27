rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_sorted_bam_dir"],
                "{sample}.sortedByCoord.bam",
            ),
            sample=config["samples"],
        ),
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_sorted_bam_dir"],
                "{sample}.sortedByCoord.bam.bai",
            ),
            sample=config["samples"],
        ),


rule sort_bam_reseq:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_aligned_bam_dir"],
            "{sample}.Aligned.bam",
        ),
    output:
        bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_sorted_bam_dir"],
            "{sample}.sortedByCoord.bam",
        ),
        index=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_sorted_bam_dir"],
            "{sample}.sortedByCoord.bam.bai",
        ),
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
