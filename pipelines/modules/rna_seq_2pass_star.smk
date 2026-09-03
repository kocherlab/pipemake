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
    log:
        "logs/star/{sample}.final.sj.log",
    params:
        index_dir=subpath(input.index_file, parent=True),
        sj_prefix=subpath(output[0], strip_suffix="SJ.out.tab"),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMmode None --outFileNamePrefix {params.sj_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} &> {log}"


rule star_pair_end_p1:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        "RNAseq/SpliceJunctions/Aligned/{sample}.SJ.out.tab",
        "RNAseq/SpliceJunctions/Aligned/{sample}.Log.final.out",
    log:
        "logs/star/{sample}.final.sj.log",
    params:
        index_dir=subpath(input.index_file, parent=True),
        sj_prefix=subpath(output[0], strip_suffix="SJ.out.tab"),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMtype None --outFileNamePrefix {params.sj_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads} &> {log}"


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
        STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --sjdbFileChrStartEnd {input.sj_file} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} &> {log}
        mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam
        """


rule star_pair_end_p2:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
        index_file="Indices/STAR/SAindex",
        sj_file="RNAseq/SpliceJunctions/Aligned/{sample}.SJ.filtered.tab",
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
        STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --sjdbFileChrStartEnd {input.sj_file} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads} &> {log}
        mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam
        """