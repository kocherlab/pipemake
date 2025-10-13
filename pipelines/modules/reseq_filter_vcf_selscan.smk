rule all:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.vcf.gz",


rule reseq_filter_vcf_bcftools:
    input:
        f"reSEQ/VCF/Unfiltered/{config['species']}_{config['assembly_version']}.vcf.gz",
    output:
        temp(
            f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.vcf.gz"
        ),
    params:
        exclude_chr=(
            f"-t ^{','.join(config['exclude_chr'])}"
            if "exclude_chr" in config and config["exclude_chr"]
            else ""
        ),
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "bcftools view -i 'F_MISSING=0.0' --min-alleles 2 --max-alleles 2 {params.exclude_chr} {input} | bcftools annotate --set-id '%CHROM\_%POS' -O z -o {output}"
