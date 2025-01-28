rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads"],
            "augustus",
            f"config.chk",
        ),


rule download_augustus_config:
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads"],
            "augustus",
            f"config.chk",
        ),
    params:
        out_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads"],
            "augustus",
        ),
        tmp_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads"],
            "tmp_augustus",
        ),
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        mkdir -p {params.tmp_dir}
        wget https://github.com/Gaius-Augustus/Augustus/archive/refs/heads/master.zip -O {params.tmp_dir}/master.zip
        unzip {params.tmp_dir}/master.zip -d {params.tmp_dir}
        mv {params.tmp_dir}/Augustus-master/config/ {params.out_dir}
        rm -rf {params.tmp_dir}
        touch {output}
        """
