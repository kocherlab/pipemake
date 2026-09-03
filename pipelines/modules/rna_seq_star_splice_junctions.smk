ruleorder: star_splice_junctions_pair_end > star_splice_junctions_single_end


rule all:
    input:
        f"RNAseq/SpliceJunctions/{config['species']}_{config['assembly_version']}.SJ.tab",


rule star_genome_generate_rnaseq:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        "Indices/STAR/SAindex",
    log:
        f"logs/star/{config['species']}_{config['assembly_version']}.index.log",
    params:
        index_dir=subpath(output[0], parent=True),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=32000,
    threads: 12
    shell:
        """
        let "index_mem_b={resources.mem_mb} * 10**6"
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir} --genomeFastaFiles {input} --limitGenomeGenerateRAM $index_mem_b --genomeSAindexNbases 13 2> {log}
        """


rule star_splice_junctions_single_end:
    input:
        r1_reads="RNAseq/FASTQs/{sample}_R1.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        temp("RNAseq/Aligned/{sample}.SJ.out.tab"),
        "RNAseq/Aligned/{sample}.Log.final.out",
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


rule star_splice_junctions_pair_end:
    input:
        r1_reads="RNAseq/FASTQs/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQs/{sample}_R2.fq.gz",
        index_file="Indices/STAR/SAindex",
    output:
        temp("RNAseq/Aligned/{sample}.SJ.out.tab"),
        "RNAseq/Aligned/{sample}.Log.final.out",
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


rule merge_splice_junctions:
    input:
        expand("RNAseq/Aligned/{sample}.SJ.out.tab", sample=config["samples"]),
    output:
        f"RNAseq/SpliceJunctions/{config['species']}_{config['assembly_version']}.SJ.tab",
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "cat {input} > {output}"
