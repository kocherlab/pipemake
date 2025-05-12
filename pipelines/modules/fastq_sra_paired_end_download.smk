rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["sra_processed_dir"],
                "{sample}_R1.fq.gz",
            ),
            sample=config["samples"],
        ),
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["sra_processed_dir"],
                "{sample}_R2.fq.gz",
            ),
            sample=config["samples"],
        ),


rule fasterq_dump_paired_end:
    output:
        r1_reads=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["sra_download_dir"],
                "{sample}_1.fastq",
            )
        ),
        r2_reads=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["sra_download_dir"],
                "{sample}_2.fastq",
            )
        ),
        tmp_dir=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["sra_download_dir"],
                "{sample}_TMP",
            )
        ),
    params:
        sra_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
        ),
        tmp_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
            "{sample}_TMP",
        ),
    singularity:
        "docker://ncbi/sra-tools:3.1.0"
    resources:
        mem_mb=16000,
        shell_exec="sh",
    threads: 4
    shell:
        """
        fasterq-dump {wildcards.sample} -O {params.sra_dir} --temp {params.tmp_dir} --threads {threads}
        rm -f {params/sra_dir}/{wildcards.sample}.fastq
        """
    
rule compress_r1_fastq:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
            "{sample}_1.fastq",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_processed_dir"],
            "{sample}_R1.fq.gz",
        ),
    params:
        sra_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
        ),
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gzip -c {input} > {output}"

rule compress_r2_fastq:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
            "{sample}_2.fastq",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_processed_dir"],
            "{sample}_R2.fq.gz",
        ),
    params:
        sra_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
        ),
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gzip -c {input} > {output}"