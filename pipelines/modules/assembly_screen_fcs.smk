rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            "fcs-adaptor",
            "fcs_adaptor_report.txt",
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
