rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_filtered_vcf_dir"],
            f"{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
        ),


rule filter_model_vcf_bcftools:
    input:
        vcf_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_unfiltered_vcf_dir"],
            f"{config['species']}_{config['assembly_version']}.vcf.gz",
        ),
        ind_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.{config['model_name']}.ind.txt",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_filtered_vcf_dir"],
            f"{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
        ),
    params:
        min_alleles=config["min_alleles"],
        max_alleles=config["max_alleles"],
        maf=config["maf_cutoff"],
        qual=config["qual_cutoff"],
        missing_cutoff=config["missing_cutoff"],
    resources:
        mem_mb=16000,
    threads: 4
    singularity:
        "docker://aewebb/bcftools:v1.20"
    shell:
        "bcftools view --samples-file {input.ind_file} {input.vcf_file} | bcftools view --min-alleles {params.min_alleles} --max-alleles {params.max_alleles} --types snps --include 'MAF>={params.maf} && QUAL>={params.qual} && F_MISSING<={params.missing_cutoff}' --output-type z --output-file {output} --threads {threads}"
