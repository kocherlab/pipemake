rule all:
    input:
        expand(
            f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca.eigenvec",
            model=config["models"],
        ),
        expand(
            f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca.eigenval",
            model=config["models"],
        ),
        expand(
            f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca.pdf",
            model=config["models"],
        ),


rule reseq_calc_freq_plink:
    input:
        bed_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bed",
        bim_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bim",
        fam_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.fam",
        ind_file=f"Models/{config['species']}.{{model}}.ind.txt",
    output:
        f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.afreq",
    params:
        bed_prefix=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned",
        pca_prefix=f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned",
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/plink2:20240418"
    shell:
        "plink2 --bfile {params.bed_prefix} --keep {input.ind_file} --freq --allow-extra-chr --out {params.pca_prefix}"


rule reseq_model_calc_pca_plink:
    input:
        bed_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bed",
        bim_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bim",
        fam_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.fam",
        ind_file=f"Models/{config['species']}.{{model}}.ind.txt",
        freq_file=f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.afreq",
    output:
        f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca.eigenvec",
        f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca.eigenval",
        f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca.log",
    params:
        bed_prefix=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned",
        pca_prefix=f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca",
        pca_count=config["pca_count"],
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/plink2:20240418"
    shell:
        "plink2 --bfile {params.bed_prefix} --keep {input.ind_file} --read-freq {input.freq_file} --pca {params.pca_count} --allow-extra-chr --out {params.pca_prefix}"


rule reseq_model_plot_pca:
    input:
        eigenvec_file=f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca.eigenvec",
        eigenval_file=f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca.eigenval",
        model_file=f"Models/{config['species']}.model",
    output:
        f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca.pdf",
    params:
        pca_dir="reSEQ/PopGen/PCA/{model}/",
        out_prefix=f"reSEQ/PopGen/PCA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.pca",
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    shell:
        "plot-pca --pca-dir {params.pca_dir} --model-file {input.model_file} --model-name {wildcards.model} --out-prefix {params.out_prefix}"
