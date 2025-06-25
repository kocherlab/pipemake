rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gtf",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
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


rule process_egapx_gtf:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "epagx",
            "epagx.gtf",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gtf",
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.no_utrs.gff",
            ),
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}",
        ),
        species=config["species"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.3"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        process-ncbi-annotations --gtf {input} --species-tag {params.species} --out-prefix {params.out_prefix}
        mv {params.out_prefix}.gff {params.out_prefix}.no_utrs.gff
        """

rule add_utrs_to_gff:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.no_utrs.gff",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        ),
    singularity:
        "docker://aewebb/ncbi-genome-tools:20250625"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "add_utrs_to_gff.py {input} > {output}"


rule gtf_to_transcripts:
    input:
        gtf_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gtf",
        ),
        assembly_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa",
        ),
    singularity:
        "docker://aewebb/gffread:v0.12.7"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gffread -x {output} -g {input.assembly_fasta} {input.gtf_file}"

rule gtf_to_proteins:
    input:
        gtf_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gtf",
        ),
        assembly_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
        ),
    singularity:
        "docker://aewebb/gffread:v0.12.7"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gffread -y {output} -g {input.assembly_fasta} {input.gtf_file}"
