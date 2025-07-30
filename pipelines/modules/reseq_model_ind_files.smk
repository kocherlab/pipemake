rule all:
    input:
        expand(f"Models/{config['species']}.{{model}}.ind.txt", model=config["models"]),
        expand(f"Models/{config['species']}.{{model}}.ind.log", model=config["models"]),


rule model_ind_file:
    input:
        f"Models/{config['species']}.model",
    output:
        f"Models/{config['species']}.{{model}}.ind.txt",
        f"Models/{config['species']}.{{model}}.ind.log",
    params:
        out_prefix=f"Models/{config['species']}.{{model}}",
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    shell:
        "model-inds --model-file {input} --model-name {wildcards.model} --out-prefix {params.out_prefix}"
