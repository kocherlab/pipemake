rule all:
    input:
        expand(
            f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.lmm.filtered.pve.txt",
            model=config["models"],
        ),
        expand(
            f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.lmm.manhattan.png",
            model=config["models"],
        ),
        expand(
            f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bslmm.param.txt",
            model=config["models"],
        ),


rule gemma_model_bed:
    input:
        bed_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bed",
        bim_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.bim",
        fam_file=f"reSEQ/PLINK/Pruned/{config['species']}_{config['assembly_version']}.pruned.fam",
        ind_file=f"Models/{config['species']}.{{model}}.ind.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bed",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bim",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.fam",
    params:
        bed_prefix=subpath(input.bed_file, strip_suffix=".bed"),
        out_prefix=subpath(output[0], strip_suffix=".bed"),
    resources:
        mem_mb=8000,
    threads: 1
    singularity:
        "docker://aewebb/plink2:20240418"
    shell:
        "plink2 --bfile {params.bed_prefix} --keep {input.ind_file} --allow-extra-chr --make-bed --out {params.out_prefix}"


rule gemma_model_phenotype_file:
    input:
        fam_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.fam",
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
        "docker://aewebb/pipemake_utils:v1.3.2"
    shell:
        "ped-phenotype-file --fam {input.fam_file} --model-file {input.model_file} --model-name {wildcards.model} --binary --out-prefix {params.out_prefix}"


rule run_gemma_gk:
    input:
        bed_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bed",
        bim_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bim",
        fam_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.fam",
        pheno_file=f"Models/GEMMA/{config['species']}.{{model}}.pheno.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.gk.cXX.txt",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.gk.log.txt",
    params:
        bed_prefix=subpath(input.bed_file, strip_suffix=".bed"),
        out_prefix=subpath(output[0], strip_suffix=".cXX.txt"),
        out_dir=subpath(output[0], parent=True),
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
        bed_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bed",
        bim_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bim",
        fam_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.fam",
        pheno_file=f"Models/GEMMA/{config['species']}.{{model}}.pheno.txt",
        gk_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.gk.cXX.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.lmm.assoc.txt",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.lmm.log.txt",
    params:
        bed_prefix=subpath(input.bed_file, strip_suffix=".bed"),
        out_prefix=subpath(output[0], basename=True, strip_suffix=".assoc.txt"),
        out_dir=subpath(output[0], parent=True),
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
        bed_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bed",
        bim_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bim",
        fam_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.fam",
        pheno_file=f"Models/GEMMA/{config['species']}.{{model}}.pheno.txt",
        gk_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.gk.cXX.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bslmm.param.txt",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bslmm.bv.txt",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bslmm.gamma.txt",
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.bslmm.hyp.txt",
    params:
        bed_prefix=subpath(input.bed_file, strip_suffix=".bed"),
        out_prefix=subpath(output[0], basename=True, strip_suffix=".param.txt"),
        out_dir=subpath(output[0], parent=True),
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
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.{{method}}.assoc.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.{{method}}.pve.txt",
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    shell:
        "calc-pve --input {input} --output {output}"


rule filter_gemma:
    input:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.lmm.pve.txt",
    output:
        filtered_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.lmm.filtered.pve.txt",
        log_file=f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.lmm.filtered.pve.txt.log",
    params:
        min_log_pvalue=config["min_log_pvalue"],
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    shell:
        "filter-gemma --gemma-file {input} --min-log-pvalue {params.min_log_pvalue} --out-filename {output.filtered_file}"


rule plot_gemma:
    input:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.lmm.assoc.txt",
    output:
        f"reSEQ/PopGen/GEMMA/{{model}}/{config['species']}_{config['assembly_version']}.lmm.manhattan.png",
    params:
        out_prefix=subpath(output[0], strip_suffix=".manhattan.png"),
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    shell:
        "manhattan-plot --input-file {input} --chrom-col chr --pos-col ps --stat-col p_wald --plot-stat-text 'Wald test p-value' --chrom-pos-sep '_' --plot-neg-log --out-prefix {params.out_prefix}"
