ruleorder: sortmerna_pair_end > sortmerna_single_end


rule all:
    input:
        "Statistics/filtered_fastq_stats.tsv",


rule sortmerna_single_end:
    input:
        r1_reads="FASTQ/Unfiltered/{sample}_R1.fq.gz",
        index_chk="Indices/sortmerna/.idx.chk",
    output:
        r1_reads="FASTQ/FIltered/rRNA/{sample}_R1.fq.gz",
        work_dir=temp(directory("FASTQ/FIltered/rRNA/{sample}")),
    params:
        read_prefix="FASTQ/FIltered/rRNA/{sample}",
        index_dir="Indices/sortmerna",
        sortmerna_db=config["sortmerna_db"],
    singularity:
        "docker://aewebb/sortmerna:v4.3.6"
    resources:
        mem_mb=48000,
    threads: 8
    shell:
        """
        sortmerna --threads {threads} --ref /opt/DBs/{params.sortmerna_db} --reads {input.r1_reads} --workdir {output.work_dir} --idx-dir {params.index_dir} --aligned {params.read_prefix} --fastx 
        mv {params.read_prefix}.fq.gz {output.r1_reads}
        """


rule sortmerna_pair_end:
    input:
        r1_reads="FASTQ/Unfiltered/{sample}_R1.fq.gz",
        r2_reads="FASTQ/Unfiltered/{sample}_R2.fq.gz",
        index_chk="Indices/sortmerna/.idx.chk",
    output:
        r1_reads="FASTQ/FIltered/rRNA/{sample}_R1.fq.gz",
        r2_reads="FASTQ/FIltered/rRNA/{sample}_R2.fq.gz",
        work_dir=temp(directory("FASTQ/FIltered/rRNA/{sample}")),
    params:
        read_prefix="FASTQ/FIltered/rRNA/{sample}",
        index_dir="Indices/sortmerna",
        sortmerna_db=config["sortmerna_db"],
    singularity:
        "docker://aewebb/sortmerna:v4.3.6"
    resources:
        mem_mb=48000,
    threads: 8
    shell:
        """
        sortmerna --threads {threads} --ref /opt/DBs/{params.sortmerna_db} --reads {input.r1_reads} --reads {input.r2_reads} --workdir {output.work_dir} --idx-dir {params.index_dir} --aligned {params.read_prefix} --fastx --out2 --paired_in
        mv {params.read_prefix}_fwd.fq.gz {output.r1_reads}
        mv {params.read_prefix}_rev.fq.gz {output.r2_reads}
        """


rule rrna_seqkit_stats:
    input:
        expand("FASTQ/FIltered/rRNA/{sample}_R1.fq.gz", sample=config["samples"]),
    output:
        "Statistics/filtered_fastq_stats.tsv",
    params:
        filtered_dir="FASTQ/FIltered/rRNA",
        filtered_wildcard="FASTQ/FIltered/rRNA/*.fq.gz",
    singularity:
        "docker://aewebb/seqkit:v2.10.0"
    resources:
        mem_mb=12000,
    threads: 8
    shell:
        "seqkit stats -j {threads} {params.filtered_wildcard} > {output} && rm -r {params.filtered_dir}"
