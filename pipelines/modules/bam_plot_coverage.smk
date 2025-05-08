rule all:
    input:
        f"Plots/Coverage/{config['species']}_{config['assembly_version']}.png",


rule plot_coverage_deeptools:
    input:
        bam=expand("Plots/BAM{sample}.sortedByCoord.bam", sample=config["samples"]),
        index=expand(
            "Plots/BAM/{sample}.sortedByCoord.bam.bai", sample=config["samples"]
        ),
    output:
        f"Plots/Coverage/{config['species']}_{config['assembly_version']}.png",
    params:
        skip_zero="--skipZeros" if config["skip_zero"] else "",
    singularity:
        "docker://aewebb/deeptools:v3.5.6"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "plotCoverage -b {input.bam} -o {output} {params.skip_zero}"
