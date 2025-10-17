rule all:
    input:
        expand(
            f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}.filtered.assoc.txt",
            model=config["models"],
            method=["lmm", "bslmm"],
        ),
        expand(
            f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}.manhattan.png",
            model=config["models"],
            method=["lmm", "bslmm"],
        ),


rule gemma_model_bed:
    input:
        bed_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bed",
        bim_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bim",
        fam_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.fam",
        ind_file=f"Models/{config['species']}.{{model}}.ind.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.bed",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.bim",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.fam",
    params:
        bed_prefix=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned",
        out_prefix=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned",
    resources:
        mem_mb=8000,
    threads: 1
    singularity:
        "docker://aewebb/plink2:20240418"
    shell:
        "plink2 --bfile {params.bed_prefix} --keep {input.ind_file} --allow-extra-chr --make-bed --out {params.out_prefix}"


rule gemma_model_phenotype_file:
    input:
        fam_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.fam",
        model_file=f"Models/{config['species']}.model",
    output:
        f"Models/GEMMA/{config['species']}.{{model}}.pheno.txt",
        f"Models/GEMMA/{config['species']}.{{model}}.pheno.log",
    params:
        out_prefix=f"Models/GEMMA/{config['species']}.{{model}}",
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    shell:
        "ped-phenotype-file --fam {input.fam_file} --model-file {input.model_file} --model-name {wildcards.model} --numeric --out-prefix {params.out_prefix}"


rule run_gemma_gk:
    input:
        bed_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.bed",
        bim_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.bim",
        fam_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.fam",
        pheno_file=f"Models/GEMMA/{config['species']}.{{model}}.pheno.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.gk.cXX.txt",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.gk.log.txt",
    params:
        bed_prefix=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned",
        out_prefix=f"{config['species']}_{config['assembly_version']}.pruned.gk",
        out_dir=f"reSEQ/PopGen/GEMMA/{{model}}",
        kinship_matrix=config["kinship_matrix"],
    resources:
        mem_mb=8000,
    threads: 1
    singularity:
        "docker://aewebb/gemma:v0.98.5"
    shell:
        "gemma -bfile {params.bed_prefix} -p {input.pheno_file} -gk {params.kinship_matrix} -outdir {params.out_dir} -o {params.out_prefix}"


rule run_gemma_lmm:
    input:
        bed_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.bed",
        bim_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.bim",
        fam_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.fam",
        pheno_file=f"Models/GEMMA/{config['species']}.{{model}}.pheno.txt",
        gk_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.gk.cXX.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.lmm.assoc.txt",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.lmm.log.txt",
    params:
        bed_prefix=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned",
        out_prefix=f"{config['species']}_{config['assembly_version']}.pruned.lmm",
        out_dir="reSEQ/PopGen/GEMMA/{model}",
        lmm_model=config["lmm_model"],
        maf=config["maf"],
    resources:
        mem_mb=8000,
    threads: 1
    singularity:
        "docker://aewebb/gemma:v0.98.5"
    shell:
        "gemma -p {input.pheno_file} -bfile {params.bed_prefix} -lmm {params.lmm_model} -k {input.gk_file} -maf {params.maf} -outdir {params.out_dir} -o {params.out_prefix}"


rule run_gemma_bslmm:
    input:
        bed_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.bed",
        bim_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.bim",
        fam_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.fam",
        pheno_file=f"Models/GEMMA/{config['species']}.{{model}}.pheno.txt",
        gk_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.gk.cXX.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.bslmm.assoc.txt",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.bslmm.log.txt",
    params:
        bed_prefix=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned",
        out_prefix=f"{config['species']}_{config['assembly_version']}.pruned.bslmm",
        out_dir="reSEQ/PopGen/GEMMA/{model}",
        bslmm_model=config["bslmm_model"],
        maf=config["maf"],
    resources:
        mem_mb=8000,
    threads: 1
    singularity:
        "docker://aewebb/gemma:v0.98.5"
    shell:
        "gemma -p {input.pheno_file} -bfile {params.bed_prefix} -bslmm {params.bslmm_model} -k {input.gk_file} -maf {params.maf} -outdir {params.out_dir} -o {params.out_prefix}"


rule calc_pve_gemma:
    input:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}.assoc.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}.pve.txt",
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.1"
    shell:
        "calc-pve --input {input} --output {output}"


rule filter_gemma:
    input:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}.pve.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}.filtered.pve.txt",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}.filtered.log",
    params:
        out_prefix=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}",
        min_log_pvalue=config["min_log_pvalue"],
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.1"
    shell:
        "filter-gemma --gemma-file {input} --min-log-pvalue {params.min_log_pvalue} --out-prefix {params.out_prefix}"


rule plot_gemma:
    input:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}.assoc.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}.manhattan.png",
    params:
        out_prefix=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.pruned.{{method}}",
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    shell:
        "manhattan-plot --input-file {input} --chrom-col chr --pos-col ps --stat-col p_wald --plot-stat-text 'Wald test p-value' --chrom-pos-sep '_' --plot-neg-log --out-prefix {params.out_prefix}"
