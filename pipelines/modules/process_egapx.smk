rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gtf",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta.gz.fai",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta.gz.gzi",
        ),


rule process_egapx_gtf:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "epagx.gtf",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gtf",
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["processed_dir"],
                f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.no_utrs.gff",
            ),
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["processed_dir"],
                f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.cvt.log",
            ),
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}",
        ),
        species=config["species"],
        assembly_version=config["assembly_version"],
        annotation_version=config["annotation_version"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.5"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        process-ncbi-annotations --gtf {input} --species-tag {params.species}_{params.assembly_version}.{params.annotation_version} --out-prefix {params.out_prefix}
        mv {params.out_prefix}.gff {params.out_prefix}.no_utrs.gff
        """


rule add_utrs_to_gff:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.no_utrs.gff",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        ),
    singularity:
        "docker://aewebb/ncbi-genome-tools:20250625"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "add_utrs_to_gff.py {input} > {output}"


rule gff_to_transcripts:
    input:
        gff_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        ),
        assembly_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["processed_dir"],
                f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa.tmp",
            )
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["processed_dir"],
                f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa.log",
            )
        ),
    singularity:
        "docker://aewebb/gffread:v0.12.7"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gffread -x {output} -g {input.assembly_fasta} {input.gff_file}"


rule update_transcripts:
    input:
        gff_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        ),
        transcript_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa.tmp",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa",
        ),
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.5"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "update-fasta --gff-file {input.gff_file} --fasta-file {input.transcript_fasta} --out-file {output}"


rule gff_to_proteins:
    input:
        gff_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        ),
        assembly_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["processed_dir"],
                f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa.tmp",
            )
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["processed_dir"],
                f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa.log",
            )
        ),
    singularity:
        "docker://aewebb/gffread:v0.12.7"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gffread -y {output} -g {input.assembly_fasta} {input.gff_file}"


rule update_proteins:
    input:
        gff_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        ),
        protein_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa.tmp",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
        ),
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.5"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "update-fasta --gff-file {input.gff_file} --fasta-file {input.protein_fasta} --out-file {output}"


rule gzip_assembly:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta.gz",
        ),
    singularity:
        "docker://aewebb/tabix:v1.11"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "bgzip -c {input} > {output}"


rule index_assembly:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta.gz",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta.gz.fai",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["processed_dir"],
            f"{config['species']}_genome_{config['assembly_version']}.fasta.gz.gzi",
        ),
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "samtools faidx {input}"
