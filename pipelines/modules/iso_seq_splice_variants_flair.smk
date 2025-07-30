ruleorder: isoseq_correct_flair_w_rnaseq > isoseq_correct_flair_w_gtf

rule all:
    input:
        os.path.join(
            "Annotations",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.splice_corrected.gtf",
        ),


rule cat_iso_seq_reads:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_fastq_dir"],
                "{sample}_R1.fq.gz",
            ),
            sample=config["samples"],
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_fastq_dir"],
                f"{config['species']}_{config['assembly_version']}.fq.gz",
            )
        ),
    shell:
        "cat {input} > {output}"

rule isoseq_align_flair:
    input:
        fasta_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
        reseq_fastq=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_fastq_dir"],
            f"{config['species']}_{config['assembly_version']}.fq.gz",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.aligned.bed",
        ),
#        fasta_file=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
#        reseq_fastq=expand("IsoSeq/FASTQ/{sample}_R1.fq.gz", sample=config["samples"]),
#    output:
#        f"Annotations/flair/{config['species']}_{config['assembly_version']}.aligned.bed",

        temp(
            f"Annotations/flair/{config['species']}_{config['assembly_version']}.aligned.bam"
        ),
        temp(
            f"Annotations/flair/{config['species']}_{config['assembly_version']}.aligned.bam.bai"
        ),
    params:
        out_prefix=f"Annotations/flair/{config['species']}_{config['assembly_version']}.aligned",
    singularity:
        "docker://aewebb/flair:v2.0.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "flair align --threads {threads} --genome {input.fasta_file} --reads {input.reseq_fastq} --output {params.out_prefix}"


rule isoseq_correct_flair_w_gtf:
    input:
        fasta_file=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
        gtf_file=f"Annotations/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.gtf",
        align_bed=f"Annotations/flair/{config['species']}_{config['assembly_version']}.aligned.bed",
    output:
        f"Annotations/flair/{config['species']}_{config['assembly_version']}_all_corrected.bed",
    params:
        out_prefix=f"Annotations/flair/{config['species']}_{config['assembly_version']}",
    singularity:
        "docker://aewebb/flair:v2.0.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "flair correct --threads {threads} --genome {input.fasta_file} --query {input.align_bed} --gtf {input.gtf_file} --output {params.out_prefix}"

rule isoseq_correct_flair_w_rnaseq:
    input:
        fasta_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
        splice_junctions=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.SJ.tab",
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
            f"{config['species']}_{config['assembly_version']}_all_corrected.bed",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}",
        ),
    singularity:
        "docker://aewebb/flair:v2.0.0"
    resources:
        mem_mb=16000,
    threads: 4
#        f"Annotations/flair/{config['species']}_{config['assembly_version']}_all_corrected.bed",
#    output:
#        directory("Annotations/flair/split_corrected"),
#    params:
#        out_prefix=f"Annotations/flair/split_corrected/{config['species']}_{config['assembly_version']}_all_corrected",
    shell:
        "flair correct --threads {threads} --genome {input.fasta_file} --query {input.align_bed} --shortread {input.splice_junctions} --output {params.out_prefix}"


rule isoseq_collapse_flair:
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
            f"{config['species']}_{config['assembly_version']}_all_corrected.bed",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.isoforms.gtf",
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                "flair",
                f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.isoforms.bed",
            ),
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["annotations_dir"],
                "flair",
                f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.isoforms.fa",
            ),
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}",
        ),
#        fasta_file=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
#        gtf_file=f"Annotations/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.gtf",
#        reseq_fastq=expand("IsoSeq/FASTQ/{sample}_R1.fq.gz", sample=config["samples"]),
#        corrected_bed=f"Annotations/flair/split_corrected/{config['species']}_{config['assembly_version']}_all_corrected.{{chrom}}.bed",
#    output:
#        "Annotations/flair/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.{{chrom}}.isoforms.gtf",
#        temp(
#            f"Annotations/flair/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.{{chrom}}.isoforms.bed"
#        ),
#        temp(
#            f"Annotations/flair/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.{{chrom}}.isoforms.fa"
#        ),
#    params:
#        out_prefix=f"Annotations/flair/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.{{chrom}}",
    singularity:
        "docker://aewebb/flair:v2.0.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "flair collapse --threads {threads} --genome {input.fasta_file} --gtf {input.gtf_file} --reads {input.reseq_fastq} --no_gtf_end_adjustment --stringent --quality -1 --check_splice --query {input.corrected_bed} --output {params.out_prefix}"


rule process_flair_output:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "flair",
            f"{config['species']}_{config['assembly_version']}.{config['annotation_version']}.isoforms.gtf",
        ),
#def aggregate_collapse(wildcards):
#    checkpoint_output = checkpoints.isoseq_split_bed_flair.get(**wildcards).output[0]
#    return expand(
#        f"Annotations/flair/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.{{chrom}}.isoforms.gtf",
#        chrom=glob_wildcards(
#            os.path.join(
#                checkpoint_output,
#                f"{config['species']}_{config['assembly_version']}_all_corrected.{{chrom}}.bed",
#            )
#        ).chrom,
#    )
#
#
#rule process_flair_output:
#    input:
#        aggregate_collapse,
    output:
        f"Annotations/{config['species']}_{config['assembly_version']}.{config['annotation_version']}.splice_corrected.gtf",
    shell:
        "cat {input} > {output}"
