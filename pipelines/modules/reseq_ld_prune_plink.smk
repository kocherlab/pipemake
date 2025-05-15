rule all:
    input:
        f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bed",
        f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bim",
        f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.fam",


rule reseq_indep_pairwise_plink:
    input:
        f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bed",
        f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bim",
        f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.fam",
    output:
        f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.prune.in",
    params:
        bed_prefix=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered",
        out_prefix=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}",
        ld_window_size=config["ld_window_size"],
        ld_window_step=config["ld_window_step"],
        ld_threshold=config["ld_threshold"],
    singularity:
        "docker://aewebb/plink2:20240418"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "plink2 --bfile {params.bed_prefix} --indep-pairwise {params.ld_window_size} {params.ld_window_step} {params.ld_threshold} --bad-ld --out {params.out_prefix} --allow-extra-chr --threads {threads}"


rule reseq_ld_prune_plink:
    input:
        bed=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bed",
        bim=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bim",
        fam=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.fam",
        prune_in=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.prune.in",
    output:
        f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bed",
        f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bim",
        f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.fam",
    params:
        input_prefix=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered",
        output_prefix=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned",
    singularity:
        "docker://aewebb/plink2:20240418"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "plink2 --bfile {params.input_prefix} --extract {input.prune_in} --make-bed --out {params.output_prefix} --allow-extra-chr --threads {threads}"
