rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff3",
        ),


rule process_braker3:
    input:
        braker3_gff=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.gff3",
        ),
        braker3_cds=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.codingseq",
        ),
        braker3_aa=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.aa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff3",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
        ),
    params:
        out_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
        ),
        species=config["species"],
        assembly_version=config["assembly_version"],
        annotation_version=config["annotation_version"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "process-braker --gff {input.braker3_gff} --fasta-cds {input.braker3_cds} --fasta-aa {input.braker3_aa} --out-dir {params.out_dir} --species {params.species} --assembly-version {params.assembly_version} --annotation-version {params.annotation_version}"
