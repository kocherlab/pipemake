rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "EggNOG",
            f"eggnog.db",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "EggNOG",
            f"eggnog.taxa.db",
        ),


rule download_eggnog_database:
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "EggNOG",
            f"eggnog.db",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "EggNOG",
            f"eggnog.taxa.db",
        ),
    params:
        data_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "EggNOG",
        ),
    singularity:
        "docker://aewebb/eggnog-mappper:v2.1.12"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "download_eggnog_data.py --data_dir {params.data_dir} -y"
