rule all:
    input:
        f"Models/{config['species']}.{config['model_name']}.ind.txt",
        f"Models/{config['species']}.{config['model_name']}.ind.log",


rule model_name_ind_file:
    input:
        f"Models/{config['species']}.model",
    output:
        f"Models/{config['species']}.{config['model_name']}.ind.txt",
        f"Models/{config['species']}.{config['model_name']}.ind.log",
    params:
        out_prefix=f"Models/{config['species']}.{config['model_name']}",
        model_name=config["model_name"],
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    shell:
        "model-inds --model-file {input} --model-name {params.model_name} --out-prefix {params.out_prefix}"
