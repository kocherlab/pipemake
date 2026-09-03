ruleorder: star_pair_end_rnaseq > star_single_end_rnaseq


rule all:
    input:
        expand("RNAseq/BAM/Aligned/{sample}.Aligned.bam", sample=config["samples"]),


rule star_single_end_rnaseq:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        "RNAseq/BAM/Aligned/{sample}.Aligned.bam",
        "RNAseq/BAM/Aligned/{sample}.Log.final.out",
    log:
        "logs/star/{sample}.aligned.log",
    params:
        index_dir=subpath(input.index_file, parent=True),
        bam_prefix=subpath(output[0], strip_suffix="Aligned.bam"),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} &> {log}
        mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam
        """


rule star_pair_end_rnaseq:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        "RNAseq/BAM/Aligned/{sample}.Aligned.bam",
        "RNAseq/BAM/Aligned/{sample}.Log.final.out",
    log:
        "logs/star/{sample}.aligned.log",
    params:
        index_dir=subpath(input.index_file, parent=True),
        bam_prefix=subpath(output[0], strip_suffix="Aligned.bam"),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads} &> {log}
        mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam
        """
