ruleorder: sortmerna_out_pair_end > sortmerna_out_single_end


rule all:
    input:
        expand("FASTQ/Filtered/{sample}_R1.fq.gz", sample=config["samples"]),


rule sortmerna_out_single_end:
    input:
        r1_reads="FASTQ/Unfiltered/{sample}_R1.fq.gz",
        index_chk="Indices/sortmerna/.idx.chk",
    output:
        r1_reads="FASTQ/Filtered/{sample}_R1.fq.gz",
        work_dir=temp(directory("FASTQ/Filtered/{sample}")),
    log:
        "logs/sortmerna/{sample}.sortmerna.log",
    params:
        read_prefix=subpath(output.r1_reads, strip_suffix="_R1.fq.gz"),
        index_dir=subpath(input.index_chk, parent=True),
        sortmerna_db=config["sortmerna_db"],
    singularity:
        "docker://aewebb/sortmerna:v4.3.7"
    resources:
        mem_mb=12000,
    threads: 8
    shell:
        """
        sortmerna --num_alignments 1 --no-best True --threads {threads} --ref /opt/DBs/{params.sortmerna_db} --reads {input.r1_reads} --workdir {output.work_dir} --idx-dir {params.index_dir} --other {params.read_prefix} --fastx &> {log}
        mv {params.read_prefix}.fq.gz {output.r1_reads}
        """


rule sortmerna_out_pair_end:
    input:
        r1_reads="FASTQ/Unfiltered/{sample}_R1.fq.gz",
        r2_reads="FASTQ/Unfiltered/{sample}_R2.fq.gz",
        index_chk="Indices/sortmerna/.idx.chk",
    output:
        r1_reads="FASTQ/Filtered/{sample}_R1.fq.gz",
        r2_reads="FASTQ/Filtered/{sample}_R2.fq.gz",
        work_dir=temp(directory("FASTQ/Filtered/{sample}")),
    log:
        "logs/sortmerna/{sample}.sortmerna.log",
    params:
        read_prefix=subpath(output.r1_reads, strip_suffix="_R1.fq.gz"),
        index_dir=subpath(input.index_chk, parent=True),
        sortmerna_db=config["sortmerna_db"],
    singularity:
        "docker://aewebb/sortmerna:v4.3.7"
    resources:
        mem_mb=12000,
    threads: 8
    shell:
        """
        sortmerna --num_alignments 1 --no-best True --threads {threads} --ref /opt/DBs/{params.sortmerna_db} --reads {input.r1_reads} --reads {input.r2_reads} --workdir {output.work_dir} --idx-dir {params.index_dir} --other {params.read_prefix} --fastx --out2 --paired_in &> {log}
        mv {params.read_prefix}_fwd.fq.gz {output.r1_reads}
        mv {params.read_prefix}_rev.fq.gz {output.r2_reads}
        """
