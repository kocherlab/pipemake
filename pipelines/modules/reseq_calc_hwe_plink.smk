rule all:
    input:
        expand(
            f"reSEQ/PopGen/HWE/{{model}}/{config['species']}_{config['assembly_version']}.hardy",
            model=config["models"],
        ),


rule reseq_calc_hwe_plink:
    input:
        bed_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bed",
        bim_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bim",
        fam_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.fam",
        ind_file=f"Models/{config['species']}.{{model}}.ind.txt",
    output:
        f"reSEQ/PopGen/HWE/{{model}}/{config['species']}_{config['assembly_version']}.hardy",
        f"reSEQ/PopGen/HWE/{{model}}/{config['species']}_{config['assembly_version']}.log",
    log:
        f"logs/plink/{config['species']}_{config['assembly_version']}.calc_hwe_{{model}}.log",
    params:
        bed_prefix=subpath(input.bed_file, strip_suffix=".bed"),
        hwe_prefix=subpath(output[0], strip_suffix=".hardy"),
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/plink2:20240418"
    shell:
        "plink2 --bfile {params.bed_prefix} --keep {input.ind_file} --hardy midp log10 --out {params.hwe_prefix} &> {log}"
