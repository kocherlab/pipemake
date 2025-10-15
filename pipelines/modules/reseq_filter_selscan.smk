rule all:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.vcf.gz",


rule reseq_filter_selscan_bcftools:
    input:
        vcf=f"reSEQ/VCF/Unfiltered/{config['species']}_{config['assembly_version']}.vcf.gz",
        ind_file=f"Models/{config['species']}.{config['model_name']}.ind.txt",
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
        "bcftools view --samples-file {input.ind_file} {params.exclude_chr} {input.vcf} | bcftools annotate --set-id '%CHROM\_%POS' -O z -o {output}"
