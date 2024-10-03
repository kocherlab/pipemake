rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.report.tsv",
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
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff3",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.proteins.fa",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.transcripts.cdna.fa",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gene-to-trans-map",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}",
        ),
    singularity:
        "docker://trinityrnaseq/trinotate:4.0.2"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        "/usr/local/src/Trinotate/util/Trinotate_GTF_or_GFF3_annot_prep.pl --annot {input.annotation_gff} --genome_fa {input.assembly_fasta) --out_prefix {params.out_prefix}"


rule create_trinotate:
    input:
        aa_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.proteins.fa",
        ),
        cds_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.transcripts.cdna.fa",
        ),
        map_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gene-to-trans-map",
        ),
    output:
        db_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.sqlite",
        ),
        data_dir=temp(
            directory(
                os.path.join(
                    config["paths"]["workflow_prefix"],
                    config["paths"]["annotations_dir"],
                    "Trinotate",
                    "Data",
                )
            )
        ),
    singularity:
        "docker://trinityrnaseq/trinotate:4.0.2"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        "/usr/local/src/Trinotate/Trinotate --create --db {output.db_file} --trinotate_data_dir {data_dir.data_dir} --use_diamond"


rule init_trinotate:
    input:
        db_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.sqlite",
        ),
        aa_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.proteins.fa",
        ),
        cds_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.transcripts.cdna.fa",
        ),
        map_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gene-to-trans-map",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                "Trinotate",
                f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.sqlite.init",
            )
        ),
    singularity:
        "docker://trinityrnaseq/trinotate:4.0.2"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        """
        /usr/local/src/Trinotate/Trinotate --db {input.db_file} -init --gene_trans_map {input.map_file} --transcript_fasta {input.cds_fasta} --transdecoder_pep {input.aa_fasta}
        touch {output}
        """


rule run_trinotate:
    input:
        db_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.sqlite",
        ),
        aa_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.proteins.fa",
        ),
        cds_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.transcripts.cdna.fa",
        ),
        init_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.sqlite.init",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                "Trinotate",
                f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.sqlite.run",
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
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.sqlite",
        ),
        run_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.sqlite.run",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "Trinotate",
            f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.report.tsv",
        ),
    singularity:
        "docker://trinityrnaseq/trinotate:4.0.2"
    resources:
        mem_mb=12000,
    threads: 1
    shell:
        "/usr/local/src/Trinotate/Trinotate --db {input.db_file} --report > {output}"
