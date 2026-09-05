ruleorder: star_pair_end_multimap > star_single_end_multimap
ruleorder: bam_to_fastq_pair_end > bam_to_fastq_single_end


rule all:
    input:
        expand(
            "RNAseq/FASTQ/Multimapped/{sample}_R1.Multimapped.fq.gz", sample=config["samples"]
        ),
        expand("RNAseq/BAM/Multimapped/{sample}.Log.final.out", sample=config["samples"]),


rule star_single_end_multimap:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        bam="RNAseq/BAM/Multimapped/{sample}.Aligned.bam",
        log="RNAseq/BAM/Multimapped/{sample}.Log.final.out",
    log:
        "logs/star/{sample}.multimapped.log",
    params:
        index_dir=subpath(input.index_file, parent=True),
        bam_prefix=subpath(output.bam, strip_suffix="Aligned.bam"),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} &> {log}
        mv {params.bam_prefix}Aligned.out.bam {output.bam}
        """


rule star_pair_end_multimap:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        bam="RNAseq/BAM/Multimapped/{sample}.Aligned.bam",
        log="RNAseq/BAM/Multimapped/{sample}.Log.final.out",
    log:
        "logs/star/{sample}.multimapped.log",
    params:
        index_dir=subpath(input.index_file, parent=True),
        bam_prefix=subpath(output.bam, strip_suffix="Aligned.bam"),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads} &> {log}
        mv {params.bam_prefix}Aligned.out.bam {output.bam}
        """


rule filter_multimapped_reads:
    input:
        "RNAseq/BAM/Multimapped/{sample}.Aligned.bam",
    output:
        "RNAseq/BAM/Multimapped/{sample}.Filtered.bam",
    log:
        "logs/samtools/{sample}.multimap_filter.log",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=4000,
    threads: 4
    shell:
        "samtools view -@ {threads} -b -e '[NH]>1' -o {output} {input} 2> {log}"


rule bam_to_fastq_single_end:
    input:
        bam="RNAseq/BAM/Multimapped/{sample}.Filtered.bam",
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
    output:
        "RNAseq/FASTQ/Multimapped/{sample}_R1.Multimapped.fq.gz",
    log:
        "logs/samtools/{sample}.multimap_fastq.log",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=4000,
    threads: 4
    shell:
        "samtools fastq -@ {threads} {input.bam} 2> {log} | gzip > {output}"


rule bam_to_fastq_pair_end:
    input:
        bam="RNAseq/BAM/Multimapped/{sample}.Filtered.bam",
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
    output:
        r1="RNAseq/FASTQ/Multimapped/{sample}_R1.Multimapped.fq.gz",
        r2="RNAseq/FASTQ/Multimapped/{sample}_R2.Multimapped.fq.gz",
    log:
        "logs/samtools/{sample}.multimap_fastq.log",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=4000,
    threads: 4
    shell:
        "samtools collate -@ {threads} -u -O {input.bam} 2> {log} | samtools fastq -@ {threads} -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n 2>> {log}"
