ruleorder: bwa_mem_pair_end_atac_seq > bwa_mem_single_end_atac_seq


rule all:
    input:
        expand("ATAC_seq/BAM/Aligned/{sample}.Aligned.bam", sample=config["samples"]),


rule bwa_index_atac_seq:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        f"Index/BWA/{config['species']}_{config['assembly_version']}.fa",
        f"Index/BWA/{config['species']}_{config['assembly_version']}.fa.amb",
        f"Index/BWA/{config['species']}_{config['assembly_version']}.fa.ann",
        f"Index/BWA/{config['species']}_{config['assembly_version']}.fa.bwt",
        f"Index/BWA/{config['species']}_{config['assembly_version']}.fa.pac",
        f"Index/BWA/{config['species']}_{config['assembly_version']}.fa.sa",
    params:
        index_fasta=f"Index/BWA/{config['species']}_{config['assembly_version']}.fa",
    singularity:
        "docker://aewebb/bwa:v0.7.18"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "cp {input} {params.index_fasta} && "
        "bwa index {params.index_fasta}"


rule bwa_mem_single_end_atac_seq:
    input:
        r1_reads="ATAC_seq/FASTQ/{sample}_R1.fastq.gz",
        index_fasta=f"Index/BWA/{config['species']}_{config['assembly_version']}.fa",
    output:
        temp("ATAC_seq/BAM/Aligned/{sample}.Aligned.bam"),
    singularity:
        "docker://aewebb/bwa:v0.7.18"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "bwa mem -t {threads} {input.index_fasta} {input.r1_reads} | samtools view --threads {threads} -bh -o {output}"


rule bwa_mem_pair_end_atac_seq:
    input:
        r1_reads="ATAC_seq/FASTQ/{sample}_R1.fastq.gz",
        r2_reads="ATAC_seq/FASTQ/{sample}_R2.fastq.gz",
        index_fasta=f"Index/BWA/{config['species']}_{config['assembly_version']}.fa",
    output:
        temp("ATAC_seq/BAM/Aligned/{sample}.Aligned.bam"),
    singularity:
        "docker://aewebb/bwa:v0.7.18"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "bwa mem -t {threads} {input.index_fasta} {input.r1_reads} {input.r2_reads} | samtools view --threads {threads} -bh -o {output}"
