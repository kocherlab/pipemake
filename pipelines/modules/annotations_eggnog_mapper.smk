rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "EggNOG",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.emapper.annotations.xlsx",
        ),


rule create_longest_aa_transcript:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "EggNOG",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "EggNOG",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep",
        ),
    singularity:
        "docker://aewebb/pipemake_utils:v1.1.1"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        "longest-transcript --input-filename {input} --output-prefix {params.out_prefix} --input-type AA --database pipemake --output-primary-id gene"


rule run_eggnog_mapper:
    input:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "EggNOG",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "EggNOG",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.emapper.annotations.xlsx",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "EggNOG",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}",
        ),
        eggnod_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "EggNOG",
        ),
        data_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "EggNOG",
        ),
    singularity:
        "docker://aewebb/eggnog-mappper:v2.1.12"
    resources:
        mem_mb=12000,
    threads: 12
    shell:
        "emapper.py --data_dir {params.data_dir} -i {input} -o {params.out_prefix} --cpu {threads} --excel --temp_dir {params.eggnod_dir}"
