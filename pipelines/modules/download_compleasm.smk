rule all:
    input:
        f"Downloads/compleasm/{config['busco_database']}.done",


rule download_compleasm_library:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        f"Downloads/compleasm/{config['busco_database']}.done",
    params:
        busco_db=config["busco_database"],
        downloads_dir="Downloads/compleasm",
    singularity:
        "docker://huangnengcsu/compleasm:v0.2.6"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "compleasm download {params.busco_db} --library_path {params.downloads_dir}"
