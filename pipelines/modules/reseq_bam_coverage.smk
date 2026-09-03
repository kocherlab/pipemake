rule all:
    input:
        expand("reSEQ/Coverage/{sample}.ln_scaled.bw", sample=config["samples"]),


rule reseq_bam_coverage_deeptools:
    input:
        bam="reSEQ/BAM/Sorted/{sample}.sortedByCoord.bam",
        index="reSEQ/BAM/Sorted/{sample}.sortedByCoord.bam.bai",
    output:
        "reSEQ/Coverage/{sample}.bw",
    log:
        "logs/deeptools/{sample}.bam_coverage.log",
    singularity:
        "docker://aewebb/deeptools:v3.5.6"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "bamCoverage -b {input.bam} -o {output} &> {log}"


rule ln_scale_coverage_wiggletools:
    input:
        "reSEQ/Coverage/{sample}.bw",
    output:
        temp("reSEQ/Coverage/{sample}.ln_scaled.wig"),
    log:
        "logs/wiggletools/{sample}.ln_scaled.log",
    singularity:
        "docker://ensemblorg/wiggletools:1.2.11"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "wiggletools ln {input} > {output} 2> {log}"


rule ln_wig_to_bigwig:
    input:
        wig_file="reSEQ/Coverage/{sample}.ln_scaled.wig",
        assembly_file=f"{config['assembly_input']}.fai",
    output:
        "reSEQ/Coverage/{sample}.ln_scaled.bw",
    log:
        "logs/wigtobigwig/{sample}.ln_scaled.log",
    singularity:
        "docker://aewebb/wigtobigwig:v2.9"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "wigToBigWig {input.wig_file} {input.assembly_file} {output} 2> {log}"
