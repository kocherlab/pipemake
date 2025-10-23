rule all:
    input:
        expand("HiFi/BAM/Sorted/{sample}.sortedByCoord.bam", sample=config["samples"]),


rule hifi_align_minimap2:
    input:
        hifi_fastq="HiFi/FASTQ/{sample}.fq.gz",
        assembly_fasta=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        bam="HiFi/BAM/Sorted/{sample}.sortedByCoord.bam",
        index="HiFi/BAM/Sorted/{sample}.sortedByCoord.bam.bai",
    singularity:
        "docker://aewebb/minimap2:v2.28"
    resources:
        mem_mb=32000,
    threads: 16
    shell:
        """
        minimap2 -ax map-hifi -t {threads} {input.assembly_fasta} {input.hifi_fastq} | samtools sort --threads {threads} -O bam -o {output.bam}
        samtools index {output.bam}
        """
