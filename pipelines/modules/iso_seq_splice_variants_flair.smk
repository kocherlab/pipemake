rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.splice_corrected.gtf",
        ),


rule isoseq_align_flair:
    input:
        fasta_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
        reseq_fastq=expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_fastq_dir"],
                "{sample}_R1.fq.gz",
            ),
            sample=config["samples"],
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.aligned.bed",
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                "flair",
                f"{config['species']}_{config['assembly_version']}.aligned.bam",
            ),
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                "flair",
                f"{config['species']}_{config['assembly_version']}.aligned.bam.bai",
            ),
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.aligned",
        ),
    singularity:
        "docker://aewebb/flair:v2.0.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "flair align --threads {threads} --genome {input.fasta_file} --reads {input.reseq_fastq} --output {params.out_prefix}"


rule isoseq_reseq_correct_flair:
    input:
        fasta_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
        gtf_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.gtf",
        ),
        align_bed=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.aligned.bed",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.corrected_all_corrected.bed",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.corrected_all_inconsistent.bed",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.corrected_cannot_verify.bed",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.corrected",
        ),
    singularity:
        "docker://aewebb/flair:v2.0.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "flair correct --threads {threads} --genome {input.fasta_file} --query {input.align_bed} --gtf {input.gtf_file} --output {params.out_prefix}"


rule isoseq_reseq_collapse_flair:
    input:
        fasta_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
        gtf_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.gtf",
        ),
        reseq_fastq=expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_fastq_dir"],
                "{sample}_R1.fq.gz",
            ),
            sample=config["samples"],
        ),
        corrected_bed=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.corrected_all_corrected.bed",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.gtf",
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                "flair",
                f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.bed",
            ),
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                "flair",
                f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.fa",
            ),
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}",
        ),
    singularity:
        "docker://aewebb/flair:v2.0.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "flair collapse --threads {threads} --genome {input.fasta_file} --gtf {input.gtf_file} --reads {input.reseq_fastq} --annotation_reliant generate --check_splice --query {input.corrected_bed} --output {params.out_prefix}"


rule process_flair_output:
    input:
        gtf_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.gtf",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.splice_corrected.gtf",
        ),
    shell:
        "cp {input.gtf_file} {output}"
