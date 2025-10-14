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


rule reseq_model_inds_pruned_plink:
    input:
        bed=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bed",
        bim=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bim",
        fam=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.fam",
        ind_file=f"Models/{config['species']}.ind.txt",
    output:
        temp(
            f"reSEQ/PLINK/GEMMA/{config['species']}_{config['assembly_version']}.pruned.bed"
        ),
        temp(
            f"reSEQ/PLINK/GEMMA/{config['species']}_{config['assembly_version']}.pruned.bim"
        ),
        temp(
            f"reSEQ/PLINK/GEMMA/{config['species']}_{config['assembly_version']}.pruned.fam"
        ),
    params:
        input_prefix=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned",
        output_prefix=f"reSEQ/PLINK/GEMMA/{config['species']}_{config['assembly_version']}.pruned",
    singularity:
        "docker://aewebb/plink2:20240418"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "plink2 --bfile {params.input_prefix} --keep {input.ind_file} --make-bed --out {params.output_prefix} --allow-extra-chr --threads {threads}"
