rule all:
    input:
        f"{config['eggnog_data_dir']}/eggnog.db",
        f"{config['eggnog_data_dir']}/eggnog.taxa.db",


rule download_eggnog_database:
    output:
        f"{config['eggnog_data_dir']}/eggnog.db",
        f"{config['eggnog_data_dir']}eggnog.taxa.db",
    params:
        data_dir=subpath(output[0], parent=True),
    log:
        "logs/eggnog-mapper/download-eggnog-database.log",
    singularity:
        "docker://aewebb/eggnog-mapper:v2.1.13"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "download_eggnog_data.py --data_dir {params.data_dir} -y &> {log}"
