rule all:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",


rule filter_basic_vcf_bcftools:
    input:
        f"reSEQ/VCF/Unfiltered/{config['species']}_{config['assembly_version']}.vcf.gz",
    output:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
    params:
        min_alleles=config["min_alleles"],
        max_alleles=config["max_alleles"],
        maf=config["maf_cutoff"],
        qual=config["qual_cutoff"],
        missing_cutoff=config["missing_cutoff"],
        include_regions=(
            f'--targets {config["include_regions"]}'
            if "include_regions" in config
            else ""
        ),
        exclude_regions=(
            f'--targets ^{config["exclude_regions"]}'
            if "exclude_regions" in config
            else ""
        ),
    resources:
        mem_mb=16000,
    threads: 4
    singularity:
        "docker://aewebb/bcftools:v1.20"
    shell:
        "bcftools view {params.include_regions} {params.exclude_regions} --min-alleles {params.min_alleles} --max-alleles {params.max_alleles} --types snps --include 'MAF>={params.maf} && QUAL>={params.qual} && F_MISSING<={params.missing_cutoff}' --output-type z --output-file {output} --threads {threads} {input}"
