rule all:
    input:
        f"Models/{config['model_name']}/{config['pop_name']}.pop",


rule pop_ind_file:
    input:
        f"Models/{config['species']}.model",
    output:
        f"Models/{config['model_name']}/{config['pop_name']}.pop",
    params:
        model_name=config["model_name"],
        out_dir=f"Models/{config['model_name']}",
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.7"
    shell:
        "model-pop-files --model-file {input} --model-name {params.model_name} --out-dir {params.out_dir}"
