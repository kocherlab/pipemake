ruleorder: sortmerna_pair_end > sortmerna_single_end


rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["filtered_fastq_dir"],
                "{sample}_R1.fq.gz",
            ),
            sample=config["samples"],
        ),


rule sortmerna_index:
    output:
        index_chk=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["index_dir"],
                "sortmerna",
                ".idx.chk",
            ),
        ),
        work_dir=temp(
            directory(
                os.path.join(
                    config["paths"]["workflow_prefix"],
                    config["paths"]["index_dir"],
                    ".sortmerna_work_dir",
                )
            ),
        ),
    params:
        index_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "sortmerna",
        ),
        sortmerna_db=config["sortmerna_db"],
    singularity:
        "docker://aewebb/sortmerna:v4.3.6"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "sortmerna -index 1 --ref /opt/DBs/{params.sortmerna_db} --idx-dir {params.index_dir} --workdir {output.work_dir} && touch {output.index_chk}"


rule sortmerna_single_end:
    input:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        index_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "sortmerna",
            ".idx.chk",
        ),
    output:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        work_dir=temp(
            directory(
                os.path.join(
                    config["paths"]["workflow_prefix"],
                    config["paths"]["filtered_fastq_dir"],
                    "{sample}",
                )
            ),
        ),
    params:
        read_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
            "{sample}",
        ),
        index_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "sortmerna",
        ),
        sortmerna_db=config["sortmerna_db"],
    singularity:
        "docker://aewebb/sortmerna:v4.3.6"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        """
        sortmerna --threads {threads} --ref /opt/DBs/{params.sortmerna_db} --reads {input.r1_reads} --workdir {output.work_dir} --idx-dir {params.index_dir} --other {params.read_prefix} --fastx 
        mv {params.read_prefix}.fq.gz {output.r1_reads}
        """


rule sortmerna_pair_end:
    input:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        r2_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}_R2.fq.gz",
        ),
        index_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "sortmerna",
            ".idx.chk",
        ),
    output:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        r2_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
            "{sample}_R2.fq.gz",
        ),
        work_dir=temp(
            directory(
                os.path.join(
                    config["paths"]["workflow_prefix"],
                    config["paths"]["filtered_fastq_dir"],
                    "{sample}",
                )
            ),
        ),
    params:
        read_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
            "{sample}",
        ),
        index_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "sortmerna",
        ),
        sortmerna_db=config["sortmerna_db"],
    singularity:
        "docker://aewebb/sortmerna:v4.3.6"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        """
        sortmerna --threads {threads} --ref /opt/DBs/{params.sortmerna_db} --reads {input.r1_reads} --reads {input.r2_reads} --workdir {output.work_dir} --idx-dir {params.index_dir} --other {params.read_prefix} --fastx --out2 --paired_in
        mv {params.read_prefix}_fwd.fq.gz {output.r1_reads}
        mv {params.read_prefix}_rev.fq.gz {output.r2_reads}
        """
