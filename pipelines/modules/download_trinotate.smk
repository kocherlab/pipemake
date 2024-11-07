rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "Trinotate",
            f"Trinotate.sqlite",
        ),


rule create_trinotate:
    output:
        db_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "Trinotate",
            f"Trinotate.sqlite",
        ),
        data_dir=directory(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["downloads_dir"],
                "Trinotate",
                "Data",
            )
        ),
    singularity:
        "docker://trinityrnaseq/trinotate:4.0.2"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        "/usr/local/src/Trinotate/Trinotate --create --db {output.db_file} --trinotate_data_dir {output.data_dir} --use_diamond"
