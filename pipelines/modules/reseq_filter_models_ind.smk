rule all:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",


rule create_models_ind_file:
    input:
        f"Models/{config['species']}.model",
    output:
        f"Models/{config['species']}.ind.txt",
    params:
        out_prefix=f"Models/{config['species']}",
        models=config["models"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.8"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "models-ind --model-file {input} --models {params.models} --out-prefix {params.out_prefix}"


rule reseq_filter_model_inds_bcftools:
    input:
        vcf=f"reSEQ/VCF/Unfiltered/{config['species']}_{config['assembly_version']}.vcf.gz",
        ind_file=f"Models/{config['species']}.ind.txt",
    output:
        temp(
            f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.vcf.gz"
        ),
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "bcftools view --samples-file {input.ind_file} {input.vcf} | bcftools view --min-alleles 2 -O z -o {output}"
