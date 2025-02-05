rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_coverage_dir"],
            f"{config['species']}_{config['assembly_version']}.png",
        ),


rule reseq_bam_coverage_deeptools:
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
            "{sample}.bw",
        ),
    singularity:
        "docker://aewebb/deeptools:v3.5.6"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "bamCoverage -b {input.bam} -o {output}"
