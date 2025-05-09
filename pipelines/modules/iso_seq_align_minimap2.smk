rule all:
    input:
        expand("IsoSeq/BAM/Aligned/{sample}.Aligned.bam", sample=config["samples"]),


rule isoseq_align_minimap2:
    input:
        fasta_file=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
        isoseq_fastqs="IsoSeq/FASTQ/{sample}_R1.fq.gz",
    output:
        temp("IsoSeq/BAM/Aligned/{sample}.Aligned.bam"),
    singularity:
        "docker://aewebb/minimap2:v2.28"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "minimap2 -t {threads} -ax splice:hq -uf {input.fasta_file} {input.isoseq_fastqs} | samtools view --threads {threads} -bh -o {output}"
