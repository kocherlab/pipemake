rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.report.tsv",
        ),


rule annotations_prep_trinotate:
    input:
        assembly_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
        annotation_gff=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.gff3",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.proteins.fa",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.transcripts.cdna.fa",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.gene-to-trans-map",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}",
        ),
    singularity:
        "docker://trinityrnaseq/trinotate:4.0.2"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        "/usr/local/src/Trinotate/util/Trinotate_GTF_or_GFF3_annot_prep.pl --annot {input.annotation_gff} --genome_fa {input.assembly_fasta} --out_prefix {params.out_prefix}"


rule copy_trinotate:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "Trinotate",
            f"Trinotate.sqlite",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.sqlite",
        ),
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "cp {input} {output}"


rule init_trinotate:
    input:
        db_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.sqlite",
        ),
        aa_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.proteins.fa",
        ),
        cds_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.transcripts.cdna.fa",
        ),
        map_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.gene-to-trans-map",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                "Trinotate",
                f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.sqlite.init",
            )
        ),
    params:
        data_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "Trinotate",
            "Data",
        ),
    singularity:
        "docker://trinityrnaseq/trinotate:4.0.2"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        """
        /usr/local/src/Trinotate/Trinotate --db {input.db_file} --init --gene_trans_map {input.map_file} --transcript_fasta {input.cds_fasta} --transdecoder_pep {input.aa_fasta} --trinotate_data_dir {params.data_dir}
        touch {output}
        """


rule run_trinotate:
    input:
        db_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.sqlite",
        ),
        aa_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.proteins.fa",
        ),
        cds_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.transcripts.cdna.fa",
        ),
        init_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.sqlite.init",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                "Trinotate",
                f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.sqlite.run",
            )
        ),
    singularity:
        "docker://trinityrnaseq/trinotate:4.0.2"
    resources:
        mem_mb=24000,
    threads: 8
    shell:
        """
        /usr/local/src/Trinotate/Trinotate --db {input.db_file} --run ALL --use_diamond --CPU {threads} --transcript_fasta {input.cds_fasta} --transdecoder_pep {input.aa_fasta}
        touch {output}
        """


rule report_trinotate:
    input:
        db_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.sqlite",
        ),
        run_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.sqlite.run",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.report.tsv",
        ),
    singularity:
        "docker://trinityrnaseq/trinotate:4.0.2"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        "/usr/local/src/Trinotate/Trinotate --db {input.db_file} --report > {output}"
