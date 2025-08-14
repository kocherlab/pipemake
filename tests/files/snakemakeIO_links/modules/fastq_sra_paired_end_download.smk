rule all:
    input:
        expand("SRA/Processed/{sample}_R1.fq.gz", sample=config["samples"]),
        expand("SRA/Processed/{sample}_R2.fq.gz", sample=config["samples"]),


rule fasterq_dump_paired_end:
    output:
        r1_reads=temp("SRA/Downloads/{sample}_1.fastq"),
        r2_reads=temp("SRA/Downloads/{sample}_2.fastq"),
    params:
        sra_dir="SRA/Downloads",
        tmp_dir="SRA/Downloads/{sample}_TMP",
    singularity:
        "docker://ncbi/sra-tools:3.1.0"
    resources:
        mem_mb=16000,
        shell_exec="sh",
    threads: 4
    shell:
        """
        fasterq-dump {wildcards.sample} -O {params.sra_dir} --temp {params.tmp_dir} --threads {threads}
        sleep 30
        rm -f {params.sra_dir}/{wildcards.sample}.fastq
        rm -rf {params.tmp_dir}
        """


rule compress_fastqs_paired_end:
    input:
        r1_reads="SRA/Downloads/{sample}_1.fastq",
        r2_reads="SRA/Downloads/{sample}_2.fastq",
    output:
        r1_reads="SRA/Processed/{sample}_R1.fq.gz",
        r2_reads="SRA/Processed/{sample}_R2.fq.gz",
    singularity:
        "docker://aewebb/pigz:v2.8"
    resources:
        mem_mb=8000,
    threads: 4
    shell:
        """
        pigz --best -c -p {threads} {input.r1_reads} > {output.r1_reads}
        pigz --best -c -p {threads} {input.r2_reads} > {output.r2_reads}
        """
