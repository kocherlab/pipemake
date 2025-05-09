rule all:
    input:
        expand("FASTQ/Unfiltered/{sample}_R1.fastq.gz", sample=config["samples"]),


rule pacbio_bam_to_fastq:
    input:
        "BAM/PacBio/{sample}.bam",
    output:
        temp("FASTQ/Unfiltered/{sample}_R1.fastq.gz"),
    singularity:
        "docker://aewebb/bamtools:v2.5.2"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "bamtools convert -format fastq -in {input} | gzip > {output}"
