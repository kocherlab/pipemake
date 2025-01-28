rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff3",
        ),


rule annotate_isoseq_braker3:
    input:
        masked_assembly=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa.masked",
        ),
        merged_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["isoseq-merged-bam-dir"],
            f"{config['species']}_{config['assembly_version']}.bam",
        ),
        protein_hints=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["homology_dir"],
            "ProteinHints.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.gff3",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.codingseq",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.aa",
        ),
    params:
        annotations_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
        ),
    singularity:
        "docker://teambraker/braker3:isoseq"
    resources:
        mem_mb=32000,
    threads: 20
    shell:
        "braker.pl --genome {input.masked_assembly} --prot_seq {input.protein_hints} --bam {input.merged_bam} -gff3 --softmasking --threads {threads} --workingdir {params.annotations_dir}"


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
        out_dir=config["paths"]["annotations_dir"],
        species=config["species"],
        assembly_version=config["assembly_version"],
        annotation_version=config["annotation_version"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.1.3"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "process-braker --gff {input.braker3_gff} --fasta-cds {input.braker3_cds} --fasta-aa {input.braker3_aa} --out-dir {params.out_dir} --species {params.species} --assembly-version {params.assembly_version} --annotation-version {params.annotation_version}"
