rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["models_dir"],
                f"{config['species']}.{{model}}.ind.txt",
            ),
            model=config["models"],
        ),
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["models_dir"],
                f"{config['species']}.{{model}}.ind.log",
            ),
            model=config["models"],
        ),


rule model_ind_file:
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
            f"{config['species']}.{{model}}.ind.txt",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.{{model}}.ind.log",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.{{model}}",
        ),
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    shell:
        "model-inds --model-file {input} --model-name {wildcards.model} --out-prefix {params.out_prefix}"
