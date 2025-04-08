rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["filtered_fastq_dir"],
                "{sample}_R1.filt.fcsfilt.fastq.gz",
            ),
            sample=config["samples"],
        ),


rule pacbio_bam_to_fasta:
    input:
        pacbio_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["pacbio_bam_dir"],
            "{sample}.bam",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["unfiltered_fastq_dir"],
                "{sample}_R1.fa",
            ),
        ),
    singularity:
        "docker://aewebb/bamtools:v2.5.2"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "bamtools convert -format fasta -in {input} -out {output}"


rule hifi_reads_screen_fcs_adaptor:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}_R1.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "fcs-adaptor",
            "fcs_adaptor_report.txt",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "fcs-adaptor",
        ),
        prok="--prok" if config["prok"] else "",
        euk="--euk" if config["euk"] else "",
    singularity:
        "docker://ncbi/fcs-adaptor:0.5.4"
    resources:
        mem_mb=16000,
        shell_exec="sh",
    threads: 1
    shell:
        "av_screen_x {input} --output {params.out_prefix} {params.prok} {params.euk}"


rule hifi_assembly_screen_hifiadapterfiltFCS:
    input:
        reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        fcs_adaptor_report=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "fcs-adaptor",
            "fcs_adaptor_report.txt",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
            "{sample}_R1.fcsfilt.fastq.gz",
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["filtered_fastq_dir"],
                "{sample}_R1.blocklist",
            )
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["filtered_fastq_dir"],
                "{sample}_R1.fastq",
            )
        ),
    params:
        out_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
        ),
    singularity:
        "docker://aewebb/hifiadapterfilt:v3.0.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "hifiadapterfiltFCS.sh -r {input.reads} -f {input.fcs_adaptor_report} -o {params.out_dir} -t {threads}"
