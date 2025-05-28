ruleorder: star_pair_end > star_single_end


rule all:
    input:
        expand(
            "RNAseq/BAM/Aligned/{sample}_R1.Unmapped.fq.gz", sample=config["samples"]
        ),
        expand("RNAseq/BAM/Aligned/{sample}.Log.final.out", sample=config["samples"]),


rule star_single_end:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        r1_reads="RNAseq/BAM/Aligned/{sample}_R1.Unmapped.fq.gz",
        log="RNAseq/BAM/Aligned/{sample}.Log.final.out",
    params:
        index_dir="Indices/STAR",
        bam_prefix="RNAseq/BAM/Aligned/{sample}.",
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outReadsUnmapped Fastx --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 1 --outSAMmode None --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads}
        gzip {params.bam_prefix}Unmapped.out.mate1
        sleep 30
        mv {params.bam_prefix}Unmapped.out.mate1.gz {output.r1_reads}
        """


rule star_pair_end:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        r1_reads="RNAseq/BAM/Aligned/{sample}_R1.Unmapped.fq.gz",
        r2_reads="RNAseq/BAM/Aligned/{sample}_R2.Unmapped.fq.gz",
        log="RNAseq/BAM/Aligned/{sample}.Log.final.out",
    params:
        index_dir="Indices/STAR",
        bam_prefix="RNAseq/BAM/Aligned/{sample}.",
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outReadsUnmapped Fastx --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads}
        gzip {params.bam_prefix}Unmapped.out.mate1
        gzip {params.bam_prefix}Unmapped.out.mate2
        sleep 30
        mv {params.bam_prefix}Unmapped.out.mate1.gz {output.r1_reads}
        mv {params.bam_prefix}Unmapped.out.mate2.gz {output.r2_reads}
        """
