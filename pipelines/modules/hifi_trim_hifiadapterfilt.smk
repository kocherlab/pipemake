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


rule hifi_reads_screen_hifiadapterfilt:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["filtered_fastq_dir"],
                "{sample}_R1.filt.fastq.gz",
            )
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
                "{sample}_R1.contaminant.blastout",
            )
        ),
    params:
        input_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}",
        ),
        output_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
        ),
        min_length=config["min_length"],
        min_match=config["min_match"],
    singularity:
        "docker://aewebb/hifiadapterfilt:v3.0.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "hifiadapterfilt.sh -p {params.input_prefix} -l {params.min_length} -m {params.min_match} -o {params.output_dir} -t {threads}"


rule hifi_assembly_screen_hifiadapterfiltFCS:
    input:
        filtered_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
            "{sample}_R1.filt.fastq.gz",
        ),
        fcs_adaptor_report=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            "fcs-adaptor",
            "fcs_adaptor_report.txt",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
            "{sample}_R1.filt.fcsfilt.fastq.gz",
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["filtered_fastq_dir"],
                "{sample}_R1.filt.blocklist",
            )
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["filtered_fastq_dir"],
                "{sample}_R1.filt.fastq",
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
        "hifiadapterfiltFCS.sh -r {input.filtered_reads} -f {input.fcs_adaptor_report} -o {params.out_dir} -t {threads}"
