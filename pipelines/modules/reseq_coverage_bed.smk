rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_coverage_dir"],
                "{sample}.bed",
            ),
            sample=config["samples"],
        ),


rule reseq_coverage_bedtools:
    input:
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
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_coverage_dir"],
            "{sample}.bed",
        ),
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "bedtools genomecov -ibam {input.bam} -bga > {output}"
