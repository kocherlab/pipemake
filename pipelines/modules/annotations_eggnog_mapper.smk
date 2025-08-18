rule all:
    input:
        f"Annotations/EggNOG/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.emapper.annotations.xlsx",


rule create_longest_aa_transcript:
    input:
        f"Annotations/{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
    output:
        f"Annotations/EggNOG/{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
    params:
        out_prefix=f"Annotations/EggNOG/{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep",
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        "longest-transcript --input-filename {input} --output-prefix {params.out_prefix} --input-type AA --database pipemake --output-primary-id gene"


rule run_eggnog_mapper:
    input:
        eggnog_db="Downloads/EggNOG/eggnog.db",
        taxa_db="Downloads/EggNOG/eggnog.taxa.db",
        transcripts=f"Annotations/EggNOG/{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
    output:
        f"Annotations/EggNOG/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.emapper.annotations.xlsx",
    params:
        out_prefix=f"Annotations/EggNOG/{config['species']}_{config['assembly_version']}.{config['annotation_version']}",
        eggnod_dir="Annotations/EggNOG/",
        data_dir="Downloads/EggNOG/",
    singularity:
        "docker://aewebb/eggnog-mappper:v2.1.12"
    resources:
        mem_mb=12000,
    threads: 12
    shell:
        "emapper.py --override --data_dir {params.data_dir} -i {input.transcripts} -o {params.out_prefix} --cpu {threads} --excel --temp_dir {params.eggnod_dir}"
