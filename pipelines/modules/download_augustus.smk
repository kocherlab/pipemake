rule all:
    input:
        f"Downloads/augustus/.config.chk",


rule download_augustus_config:
    output:
        "Downloads/augustus/.config.chk",
    params:
        out_dir="Downloads/augustus",
        tmp_dir="Downloads/tmp_augustus",
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
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
