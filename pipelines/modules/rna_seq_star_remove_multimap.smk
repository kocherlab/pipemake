ruleorder: star_pair_end_rnaseq > star_single_end_rnaseq


rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["rnaseq_aligned_bam_dir"],
                "{sample}.Aligned.bam",
            ),
            sample=config["samples"],
        ),


rule star_genome_generate_rnaseq:
    input:
        fasta_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
    output:
        index_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "STAR",
            "SAindex",
        ),
    params:
        index_dir=directory(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["index_dir"],
                "STAR",
            )
        ),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=32000,
    threads: 12
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir} --genomeFastaFiles {input.fasta_file} --limitGenomeGenerateRAM $index_mem_b --genomeSAindexNbases 13"


rule star_single_end_rnaseq:
    input:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        index_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "STAR",
            "SAindex",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_aligned_bam_dir"],
            "{sample}.Aligned.bam",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_aligned_bam_dir"],
            "{sample}.Log.final.out",
        ),
    params:
        index_dir=os.path.join(
            config["paths"]["workflow_prefix"], config["paths"]["index_dir"], "STAR"
        ),
        bam_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_aligned_bam_dir"],
            "{sample}.",
        ),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    resources:
        mem_mb=1,
    shell:
        """
        STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads}
        mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam
        """


rule star_pair_end_rnaseq:
    input:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        r2_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_fastq_dir"],
            "{sample}_R2.fq.gz",
        ),
        index_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "STAR",
            "SAindex",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_aligned_bam_dir"],
            "{sample}.Aligned.bam",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_aligned_bam_dir"],
            "{sample}.Log.final.out",
        ),
    params:
        index_dir=os.path.join(
            config["paths"]["workflow_prefix"], config["paths"]["index_dir"], "STAR"
        ),
        bam_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_aligned_bam_dir"],
            "{sample}.",
        ),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads}
        mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam
        """
