rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["unfiltered_fastq_dir"],
                "{sample}_R1.fq.gz",
            ),
            sample=config["samples"],
        ),


rule pacbio_bam_to_fastq:
    input:
        pacbio_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["pacbio_bam_dir"],
            "{sample}.bam",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["unfiltered_fastq_dir"],
                "{sample}_R1.fq.gz",
            ),
        ),
    singularity:
        "docker://aewebb/bamtools:v2.5.2"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "bamtools convert -format fastq -in {input}| gzip > {output}"
