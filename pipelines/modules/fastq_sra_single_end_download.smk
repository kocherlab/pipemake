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


rule fasterq_dump_single_end:
    output:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
            "{sample}_R1.fq.gz",
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
        if [ -f {params.sra_dir}/{wildcards.sample}_2.fastq ]; then
            echo "Paired-end reads detected for {wildcards.sample}. Exiting."
            exit 1
        fi
        rm {params.sra_dir}/{wildcards.sample}.fastq
        rm -r {params.tmp_dir}
        gzip {params.sra_dir}/*.fastq
        mv {params.sra_dir}/{wildcards.sample}_1.fastq.gz {output.r1_reads}
        """


rule fastq_sra_process_reads:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_download_dir"],
            "{SRA}_R1.fq.gz",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["sra_processed_dir"],
            "{SRA}_R1.fq.gz",
        ),
    resources:
        mem_mb=1000,
    threads: 1
    run:
        import os

        os.symlink(os.path.abspath(input), os.path.abspath(output))
