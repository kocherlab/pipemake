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
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
            "{sample}_R1.fq.gz",
        ),
        r2_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
            "{sample}_R2.fq.gz",
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
        mem_mb=8000,
        shell_exec="sh",
    threads: 1
    shell:
        """
        fasterq-dump {wildcards.sample} -O {params.sra_dir} --temp {params.tmp_dir}
        rm {params.sra_dir}/{wildcards.sample}.fastq
        rm -r {params.tmp_dir}
        gzip {params.sra_dir}/*.fastq
        mv {params.sra_dir}/{wildcards.sample}_1.fastq.gz {output.r1_reads}
        mv {params.sra_dir}/{wildcards.sample}_2.fastq.gz {output.r2_reads}
        """


rule fastq_sra_process_reads:
    input:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
            "{SRA}_R1.fq.gz",
        ),
        r2_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
            "{SRA}_R2.fq.gz",
        ),
    output:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_processed_dir"],
            "{SRA}_R1.fq.gz",
        ),
        r2_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_processed_dir"],
            "{SRA}_R2.fq.gz",
        ),
    params:
        sra_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
            "{SRA}_DIR",
        ),
    resources:
        mem_mb=1000,
    threads: 1
    run:
        import os

        os.symlink(os.path.abspath(input.r1_reads), os.path.abspath(output.r1_reads))
        os.symlink(os.path.abspath(input.r2_reads), os.path.abspath(output.r2_reads))
