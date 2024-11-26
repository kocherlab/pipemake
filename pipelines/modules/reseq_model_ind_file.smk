rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.{config['model_name']}.ind.txt",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.{config['model_name']}.ind.log",
        ),


rule model_name_ind_file:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.model",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.{config['model_name']}.ind.txt",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.{config['model_name']}.ind.log",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.{config['model_name']}",
        ),
        model_name=config["model_name"],
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.1.3"
    shell:
        "model-inds --model-file {input} --model-name {params.model_name} --out-prefix {params.out_prefix}"
