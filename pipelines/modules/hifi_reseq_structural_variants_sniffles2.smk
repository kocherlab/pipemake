rule all:
    input:
        expand("reSEQ/VCFs/{sample}.vcf.gz", sample=config["samples"]),
        expand("reSEQ/VCFs/{sample}.snf", sample=config["samples"]),


rule hifi_structural_variants_sniffles2:
    input:
        hifi_bam="HiFi/BAM/Sorted/{sample}.sortedByCoord.bam",
        assembly_fasta=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        vcf="reSEQ/VCFs/{sample}.vcf.gz",
        sif="reSEQ/VCFs/{sample}.snf",
    params:
        tandem_repeats=(
            f"--tandem-repeats {config['tandem_repeats']}"
            if config.get("tandem_repeats")
            else ""
        ),
    singularity:
        "docker://aewebb/sniffles:v2.6.3"
    resources:
        mem_mb=4000,
    threads: 8
    shell:
        "sniffles --input {input.hifi_bam} --reference {input.assembly_fasta} --vcf {output.vcf} --snf {output.sif} --threads {threads} {params.tandem_repeats}"
