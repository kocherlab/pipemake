rule all:
    input:
        expand("reSEQ/Coverage/{sample}.ln_scaled.bw", sample=config["samples"]),


rule reseq_bam_coverage_deeptools:
    input:
        bam="reSEQ/BAM/Sorted/{sample}.sortedByCoord.bam",
        index="reSEQ/BAM/Sorted/{sample}.sortedByCoord.bam.bai",
    output:
        "reSEQ/Coverage/{sample}.bw",
    singularity:
        "docker://aewebb/deeptools:v3.5.6"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "bamCoverage -b {input.bam} -o {output}"


rule ln_scale_coverage_wiggletools:
    input:
        "reSEQ/Coverage/{sample}.bw",
    output:
        temp("reSEQ/Coverage/{sample}.ln_scaled.wig"),
    singularity:
        "docker://ensemblorg/wiggletools:1.2.11"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "wiggletools ln {input} > {output}"


rule ln_wig_to_bigwig:
    input:
        wig_file="reSEQ/Coverage/{sample}.ln_scaled.wig",
        assembly_file=f"Assembly/{config['species']}_{config['assembly_version']}.fa.fai",
    output:
        "reSEQ/Coverage/{sample}.ln_scaled.bw",
    singularity:
        "docker://aewebb/wigtobigwig:v2.9"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "wigToBigWig {input.wig_file} {input.assembly_file} {output}"
