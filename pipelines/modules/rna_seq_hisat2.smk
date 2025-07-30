ruleorder: hisat2_pair_end > hisat2_single_end


rule all:
    input:
        expand("RNAseq/BAM/Aligned/{sample}.Aligned.bam", sample=config["samples"]),


rule hisat2_build:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        "Indices/hisat2/hisat2_index.1.ht2",
        "Indices/hisat2/hisat2_index.2.ht2",
        "Indices/hisat2/hisat2_index.3.ht2",
        "Indices/hisat2/hisat2_index.4.ht2",
        "Indices/hisat2/hisat2_index.5.ht2",
        "Indices/hisat2/hisat2_index.6.ht2",
        "Indices/hisat2/hisat2_index.7.ht2",
        "Indices/hisat2/hisat2_index.8.ht2",
    params:
        index_prefix="Indices/hisat2/hisat2_index",
    singularity:
        "docker://aewebb/hisat2:v2.2.1"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "hisat2-build -p {threads} {input} {params.index_prefix}"


rule hisat2_single_end:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        index1="Indices/hisat2/hisat2_index.1.ht2",
    output:
        "RNAseq/BAM/Aligned/{sample}.Aligned.bam",
    params:
        index_prefix="Indices/hisat2/hisat2_index",
    singularity:
        "docker://aewebb/hisat2:v2.2.1"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "hisat2 --threads {threads} --dta -q -x {params.index_prefix} -U {input.r1_reads} | samtools view -@ {threads} -bh -o {output}"


rule hisat2_pair_end:
    input:
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
        index1="Indices/hisat2/hisat2_index.1.ht2",
    output:
        "RNAseq/BAM/Aligned/{sample}.Aligned.bam",
    params:
        index_prefix="Indices/hisat2/hisat2_index",
    singularity:
        "docker://aewebb/hisat2:v2.2.1"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "hisat2 --threads {threads} --dta -q -x {params.index_prefix} -1 {input.r1_reads} -2 {input.r2_reads} | samtools view -@ {threads} -bh -o {output}"
