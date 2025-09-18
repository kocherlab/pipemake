rule all:
    input:
        expand(
            f"reSEQ/PopGen/{config['species']}.{{model}}.report",
            model=config["models"],
        ),


rule reseq_model_fst_phenotype_file:
    input:
        fam_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.fam",
        model_file=f"Models/{config['species']}.model",
    output:
        f"Models/Fst/{config['species']}.{{model}}.pheno.txt",
        f"Models/Fst/{config['species']}.{{model}}.pheno.log",
    params:
        out_prefix=f"Models/Fst/{config['species']}.{{model}}",
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    shell:
        "ped-phenotype-file --fam {input.fam_file} --model-file {input.model_file} --model-name {wildcards.model} --out-prefix {params.out_prefix} --out-format plink2 --pheno-header {wildcards.model}"


checkpoint reseq_model_calc_fst_plink:
    input:
        bed_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bed",
        bim_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bim",
        fam_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.fam",
        pheno_file=f"Models/Fst/{config['species']}.{{model}}.pheno.txt",
        ind_file=f"Models/{config['species']}.{{model}}.ind.txt",
    output:
        fst_dir=directory("reSEQ/PopGen/Fst/{model}"),
    params:
        bed_prefix=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered",
        fst_prefix=f"reSEQ/PopGen/Fst/{{model}}/{config['species']}_{config['assembly_version']}.filtered",
        fst_method=config["fst_method"],
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/plink2:20240418"
    shell:
        """
        mkdir -p {output.fst_dir}
        plink2 --bfile {params.bed_prefix} --pheno {input.pheno_file} --keep {input.ind_file} --fst {wildcards.model} report-variants method={params.fst_method} --allow-extra-chr --out {params.fst_prefix}
        """


def get_fst_files(wildcards):
    checkpoint_output = checkpoints.reseq_model_calc_fst_plink.get(**wildcards).output[
        "fst_dir"
    ]
    return expand(
        os.path.join(
            checkpoint_output,
            f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst.var",
        ),
        pair=glob_wildcards(
            os.path.join(
                checkpoint_output,
                f"{config['species']}_{config['assembly_version']}.filtered.{{pair}}.fst.var",
            )
        ).pair,
    )


rule reseq_model_fst_tmp_report:
    input:
        get_fst_files,
    output:
        temp(f"reSEQ/PopGen/{config['species']}.{{model}}.report"),
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "echo {input} > {output}"
