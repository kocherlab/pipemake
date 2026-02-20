rule all:
    input:
        f"Annotations/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff3",


rule process_braker3:
    input:
        braker3_gff="Annotations/BRAKER3/braker.gff3",
        braker3_cds="Annotations/BRAKER3/braker.codingseq",
        braker3_aa="Annotations/BRAKER3/braker.aa",
    output:
        f"Annotations/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff3",
        f"Annotations/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa",
        f"Annotations/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
    params:
        out_dir="Annotations",
        species=config["species"],
        assembly_version=config["assembly_version"],
        annotation_version=config["annotation_version"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.7"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "process-braker --gff {input.braker3_gff} --fasta-cds {input.braker3_cds} --fasta-aa {input.braker3_aa} --out-dir {params.out_dir} --species {params.species} --assembly-version {params.assembly_version} --annotation-version {params.annotation_version}"
