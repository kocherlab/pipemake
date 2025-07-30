ruleorder: star_pair_end_p1 > star_single_end_p1
ruleorder: star_pair_end_p2 > star_single_end_p2


rule all:
    input:
        expand("RNAseq/BAM/Aligned/{sample}.Aligned.bam", sample=config["samples"]),


rule star_single_end_p1:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        "RNAseq/SpliceJunctions/Aligned/{sample}.SJ.out.tab",
        "RNAseq/SpliceJunctions/Aligned/{sample}.Log.final.out",
    params:
        index_dir="Indices/STAR",
        sj_prefix="RNAseq/SpliceJunctions/Aligned/{sample}.",
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMmode None --outFileNamePrefix {params.sj_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads}"


rule star_pair_end_p1:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        "RNAseq/SpliceJunctions/Aligned/{sample}.SJ.out.tab",
        "RNAseq/SpliceJunctions/Aligned/{sample}.Log.final.out",
    params:
        index_dir="Indices/STAR",
        sj_prefix="RNAseq/SpliceJunctions/Aligned/{sample}.",
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMtype None --outFileNamePrefix {params.sj_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads}"


rule filter_star_sj_file:
    input:
        "RNAseq/SpliceJunctions/Aligned/{sample}.SJ.out.tab",
    output:
        "RNAseq/SpliceJunctions/Aligned/{sample}.SJ.filtered.tab",
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "cat {input} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output}"


rule star_single_end_p2:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        index_file="Indices/STAR/SAindex",
        sj_file="RNAseq/SpliceJunctions/Aligned/{sample}.SJ.filtered.tab",
    output:
        temp("RNAseq/BAM/Aligned/{sample}.Aligned.bam"),
        "RNAseq/BAM/Aligned/{sample}.Log.final.out",
    params:
        index_dir="Indices/STAR",
        bam_prefix="RNAseq/BAM/Aligned/{sample}.",
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    resources:
        mem_mb=1,
    shell:
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --sjdbFileChrStartEnd {input.sj_file} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} && "
        "mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam"


rule star_pair_end_p2:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
        index_file="Indices/STAR/SAindex",
        sj_file="RNAseq/SpliceJunctions/Aligned/{sample}.SJ.filtered.tab",
    output:
        "RNAseq/BAM/Aligned/{sample}.Aligned.bam",
        "RNAseq/BAM/Aligned/{sample}.Log.final.out",
    params:
        index_dir="Indices/STAR",
        bam_prefix="RNAseq/BAM/Aligned/{sample}.",
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --sjdbFileChrStartEnd {input.sj_file} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads} && "
        "mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam"
