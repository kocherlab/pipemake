rule all:
    input:
        f"Annotations/EggNOG/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.emapper.annotations.xlsx",


rule create_longest_aa_transcripts:
    input:
        f"Annotations/Unprocessed/{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
    output:
        f"Annotations/EggNOG/{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
    log:
        f"logs/create-longest-aa-transcripts/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.pep.log",
    params:
        out_prefix=subpath(output[0], strip_suffix=".fa"),
    singularity:
        "docker://aewebb/pipemake_utils:v1.5.0"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        "longest-transcript --input-filename {input} --output-prefix {params.out_prefix} --input-type AA --database pipemake --output-primary-id gene &> {log}"


rule run_eggnog_mapper:
    input:
        eggnog_db=f"{config['eggnog_data_dir']}/eggnog.db",
        taxa_db=f"{config['eggnog_data_dir']}/eggnog.taxa.db",
        transcripts=f"Annotations/EggNOG/{config['species']}_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
    output:
        f"Annotations/EggNOG/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.emapper.annotations.xlsx",
    log:
        f"logs/eggnog-mapper/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.log",
    params:
        out_prefix=subpath(output[0], strip_suffix=".emapper.annotations.xlsx"),
        eggnod_dir=subpath(output[0], parent=True),
        data_dir=subpath(input.eggnog_db, parent=True),
    singularity:
        "docker://aewebb/eggnog-mapper:v2.1.13"
    resources:
        mem_mb=12000,
    threads: 12
    shell:
        "emapper.py --override --data_dir {params.data_dir} -i {input.transcripts} -o {params.out_prefix} --cpu {threads} --excel --temp_dir {params.eggnod_dir} &> {log}"
