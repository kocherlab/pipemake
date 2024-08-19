rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gff3",
        ),


rule annotate_3prime_utrs_peaks2utr:
    input:
        merged_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_bam_dir"],
            f"{config['species']}_{config['assembly_version']}.bam",
        ),
        base_gff=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.gff3",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "peaks2utr",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gff3",
        ),
    resources:
        mem_mb=120000,
    threads: 10
    singularity:
        "docker://aewebb/peaks2utr:v1.2.5"
    shell:
        "peaks2utr {input.base_gff} {input.merged_bam} -o {output} -p {threads}"


rule process_peaks2utr_utrs:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "peaks2utr",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gff3",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gff3",
        ),
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "grep -v '##sequence-region' {input} | grep -v '###' > {output}"
