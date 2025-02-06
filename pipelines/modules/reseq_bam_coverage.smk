rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_coverage_dir"],
                "{sample}.ln_scaled.bw",
            ),
            sample=config["samples"],
        ),


rule reseq_bam_coverage_deeptools:
    input:
        bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_sorted_bam_dir"],
            "{sample}.sortedByCoord.bam",
        ),
        index=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_sorted_bam_dir"],
            "{sample}.sortedByCoord.bam.bai",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_coverage_dir"],
            "{sample}.bw",
        ),
    singularity:
        "docker://aewebb/deeptools:v3.5.6"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "bamCoverage -b {input.bam} -o {output}"


rule ln_scale_coverage_wiggletools:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_coverage_dir"],
            "{sample}.bw",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_coverage_dir"],
                "{sample}.ln_scaled.wig",
            ),
        ),
    singularity:
        "docker://ensemblorg/wiggletools:1.2.11"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "wiggletools ln {input} > {output}"


rule ln_wig_to_bigwig:
    input:
        wig_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_coverage_dir"],
            "{sample}.ln_scaled.wig",
        ),
        assembly_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa.fai",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_coverage_dir"],
            "{sample}.ln_scaled.bw",
        ),
    singularity:
        "docker://aewebb/wigtobigwig:v2.9"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "wigToBigWig {input.wig_file} {input.assembly_file} {output}"
