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


rule assembly_screen_fcs_adaptor:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            "fcs-adaptor",
            "fcs_adaptor_report.txt",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            "fcs-adaptor",
        ),
        prok="--prok" if "prok" in config else "",
        euk="--euk" if "euk" in config else "",
    singularity:
        "docker://ncbi/fcs-adaptor:0.5.4"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "/app/fcs/bin/av_screen_x {input} --output {params.out_prefix} {params.prok} {params.euk}"


rule hifi_screen_hifiadapterfilt:
    input:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        fcs_adaptor_report=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            "fcs-adaptor",
            "fcs_adaptor_report.txt",
        ),
    output:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    params:
        output_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
        ),
    singularity:
        "docker://aewebb/fastp:v0.23.4"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "hifiadapterfiltFCS.sh -r {input.r1_reads} -f {input.fcs_adaptor_report} -o {params.output_prefix} -t {threads}"
