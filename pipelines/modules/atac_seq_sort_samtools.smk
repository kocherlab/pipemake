rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["atacseq_sorted_bam_dir"],
                "{sample}.sortedByCoord.bam",
            ),
            sample=config["samples"],
        ),


rule sort_bam:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_aligned_bam_dir"],
            "{sample}.Aligned.bam",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_sorted_bam_dir"],
            "{sample}.sortedByCoord.bam",
        ),
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output} {input}"
