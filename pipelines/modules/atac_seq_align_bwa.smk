ruleorder: bwa_mem_pair_end_atac_seq > bwa_mem_single_end_atac_seq


rule all:
    input:
        expand("ATAC_seq/BAM/Aligned/{sample}.Aligned.bam", sample=config["samples"]),


rule bwa_index_atac_seq:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        fasta=f"Index/BWA/{config['species']}_{config['assembly_version']}.fa",
        index=f"Index/BWA/{config['species']}_{config['assembly_version']}.fa.bwt.2bit.64",
    log:
        f"logs/bwa-mem2/{config['species']}_{config['assembly_version']}.index.log",
    singularity:
        "docker://aewebb/bwa-mem2:v2.2.1"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "cp {input} {output.fasta} && bwa-mem2 index {output.fasta} 2> {log}"


rule bwa_mem_single_end_atac_seq:
    input:
        r1_reads="ATAC_seq/FASTQ/{sample}_R1.fastq.gz",
        index_fasta=f"Index/BWA/{config['species']}_{config['assembly_version']}.fa",
    output:
        temp("ATAC_seq/BAM/Aligned/{sample}.Aligned.bam"),
    log:
        "logs/bwa-mem2/{sample}.aligned.log",
    singularity:
        "docker://aewebb/bwa-mem2:v2.2.1"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "bwa-mem2 mem -t {threads} {input.index_fasta} {input.r1_reads} | samtools view --threads {threads} -bh -o {output} 2> {log}"


rule bwa_mem_pair_end_atac_seq:
    input:
        r1_reads="ATAC_seq/FASTQ/{sample}_R1.fastq.gz",
        r2_reads="ATAC_seq/FASTQ/{sample}_R2.fastq.gz",
        index_fasta=f"Index/BWA/{config['species']}_{config['assembly_version']}.fa",
    output:
        temp("ATAC_seq/BAM/Aligned/{sample}.Aligned.bam"),
    log:
        "logs/bwa-mem2/{sample}.aligned.log",
    singularity:
        "docker://aewebb/bwa-mem2:v2.2.1"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "bwa-mem2 mem -t {threads} {input.index_fasta} {input.r1_reads} {input.r2_reads} | samtools view --threads {threads} -bh -o {output} 2> {log}"
