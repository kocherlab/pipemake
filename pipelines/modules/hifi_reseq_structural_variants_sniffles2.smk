rule all:
    input:
        expand("reSEQ/VCFs/{sample}.SVs.vcf.gz", sample=config["samples"]),
        f"reSEQ/VCFs/{config['species']}_{config['assembly_version']}.SVs.vcf.gz",


rule hifi_structural_variants_sniffles2:
    input:
        hifi_bam="HiFi/BAM/Sorted/{sample}.sortedByCoord.bam",
        assembly_fasta=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        vcf="reSEQ/VCFs/{sample}.SVs.vcf.gz",
        sif=temp("reSEQ/SNFs/{sample}.snf"),
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


rule multisample_vcf_sniffles2:
    input:
        expand("reSEQ/SNFs/{sample}.snf", sample=config["samples"]),
    output:
        f"reSEQ/VCFs/{config['species']}_{config['assembly_version']}.SVs.vcf.gz",
    singularity:
        "docker://aewebb/sniffles:v2.6.3"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "sniffles --input {input} --vcf {output}"
