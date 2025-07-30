ruleorder: star_pair_end_rnaseq > star_single_end_rnaseq


rule all:
    input:
        expand("RNAseq/BAM/Aligned/{sample}.Aligned.bam", sample=config["samples"]),


rule star_genome_generate_rnaseq:
    input:
        fasta_file=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        index_file="Indices/STAR/SAindex",
    params:
        index_dir="Indices/STAR",
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=32000,
    threads: 12
    shell:
        """
        let "index_mem_b={resources.mem_mb} * 10**6"
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir} --genomeFastaFiles {input.fasta_file} --limitGenomeGenerateRAM $index_mem_b --genomeSAindexNbases 13
        """


rule star_single_end_rnaseq:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        index_file="Indices/STAR/SAindex",
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
    resources:
        mem_mb=1,
    shell:
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} && "
        "mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam"


rule star_pair_end_rnaseq:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
        index_file="Indices/STAR/SAindex",
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
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads} && "
        "mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam"
