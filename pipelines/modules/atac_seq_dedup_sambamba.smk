rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["atacseq_dedup_bam_dir"],
                "{sample}.deduplicated.bam",
            ),
            sample=config["samples"],
        ),


rule dedup_bam_sambamba:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_sorted_bam_dir"],
            "{sample}.sortedByCoord.bam",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_dedup_bam_dir"],
            "{sample}.deduplicated.bam",
        ),
    singularity:
        "docker://aewebb/sambamba:v1.0.1"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "sambamba markdup -t {threads} {input} {output}"
