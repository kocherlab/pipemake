rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_coverage_dir"],
            f"{config['species']}_{config['assembly_version']}.png",
        ),


rule reseq_coverage_deeptools:
    input:
        bam=expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_sorted_bam_dir"],
                "{sample}.sortedByCoord.bam",
            ),
            sample=config["samples"],
        ),
        index=expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_sorted_bam_dir"],
                "{sample}.sortedByCoord.bam.bai",
            ),
            sample=config["samples"],
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_coverage_dir"],
            f"{config['species']}_{config['assembly_version']}.png",
        ),
    singularity:
        "docker://aewebb/deeptools:v3.5.6"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "plotCoverage -b {input.bam} -o {output}"
