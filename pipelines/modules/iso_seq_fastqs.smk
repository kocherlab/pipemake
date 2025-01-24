rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["isoseq_fastq_dir"],
                "{sample}.fastq.gz",
            ),
            sample=config["samples"],
        ),


rule pacbio_to_fastq:
    input:
        pacbio_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["pacbio_bam_dir"],
            "{sample}.bam",
        ),
        pacbio_idx=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["pacbio_bam_dir"],
            "{sample}.bam.pbi",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["isoseq_fastq_dir"],
                "{sample}.fastq.gz",
            ),
        ),
    params:
        pacbio_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["pacbio_bam_dir"],
            "{sample}",
        ),
    singularity:
        "docker://aewebb/pbtk:v3.4.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "bam2fastq --num-threads {threads} --output {params.pacbio_prefix} {input.pacbio_bam}"
