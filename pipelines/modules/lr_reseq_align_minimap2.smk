rule all:
    input:
        expand("LR_reSEQ/BAM/Aligned/{sample}.Aligned.bam", sample=config["samples"]),


rule longread_reseq_align_minimap2:
    input:
        fasta_file=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
        reseq_fastq="LR_reSEQ/FASTQ/{sample}.fq.gz",
    output:
        temp("LR_reSEQ/BAM/Aligned/{sample}.Aligned.bam"),
    singularity:
        "docker://aewebb/minimap2:v2.28"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "minimap2 -t {threads} -ax map-hifi -uf {input.fasta_file} {input.reseq_fastq} | samtools view --threads {threads} -bh -o {output}"
