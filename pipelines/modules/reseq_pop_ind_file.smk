rule all:
    input:
        f"Models/{config['model_name']}/{config['pop_name']}.pop",


rule pop_ind_file:
    input:
        f"Models/{config['species']}.model",
    output:
        f"Models/{config['model_name']}/{config['pop_name']}.pop",
    log:
        f"logs/model-pop-files/{config['model_name']}_{config['pop_name']}.log",
    params:
        model_name=config["model_name"],
        out_dir=subpath(output[0], parent=True),
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.4.3"
    shell:
        "model-pop-files --model-file {input} --model-name {params.model_name} --out-dir {params.out_dir} &> {log}"
