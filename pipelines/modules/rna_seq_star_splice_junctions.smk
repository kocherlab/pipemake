ruleorder: star_splice_junctions_pair_end > star_splice_junctions_single_end


rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_splice_aligned_dir"],
            f"{config['species']}_{config['assembly_version']}.SJ.tab",
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
        """
        let "index_mem_b={resources.mem_mb} * 10**6"
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir} --genomeFastaFiles {input.fasta_file} --limitGenomeGenerateRAM $index_mem_b --genomeSAindexNbases 13
        """

rule star_splice_junctions_single_end:
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
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["rnaseq_splice_aligned_dir"],
                "{sample}.SJ.out.tab",
            )
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_splice_aligned_dir"],
            "{sample}.Log.final.out",
        ),
    params:
        index_dir=os.path.join(
            config["paths"]["workflow_prefix"], config["paths"]["index_dir"], "STAR"
        ),
        sj_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_splice_aligned_dir"],
            "{sample}.",
        ),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMmode None --outFileNamePrefix {params.sj_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads}"


rule star_splice_junctions_pair_end:
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
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["rnaseq_splice_aligned_dir"],
                "{sample}.SJ.out.tab",
            )
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_splice_aligned_dir"],
            "{sample}.Log.final.out",
        ),
    params:
        index_dir=os.path.join(
            config["paths"]["workflow_prefix"], config["paths"]["index_dir"], "STAR"
        ),
        sj_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_splice_aligned_dir"],
            "{sample}.",
        ),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMtype None --outFileNamePrefix {params.sj_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads}"

rule merge_splice_junctions:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["rnaseq_splice_aligned_dir"],
                "{sample}.SJ.out.tab",
            ),
            sample=config["samples"],
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_splice_aligned_dir"],
            f"{config['species']}_{config['assembly_version']}.SJ.tab",
        ),
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "cat {input} > {output}"