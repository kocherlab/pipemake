rule all:
    input:
        "Downloads/EggNOG/eggnog.db",
        "Downloads/EggNOG/eggnog.taxa.db",


rule download_eggnog_database:
    output:
        "Downloads/EggNOG/eggnog.db",
        "Downloads/EggNOG/eggnog.taxa.db",
    params:
        data_dir="Downloads/EggNOG",
    singularity:
        "docker://aewebb/eggnog-mappper:v2.1.12"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "download_eggnog_data.py --data_dir {params.data_dir} -y"
