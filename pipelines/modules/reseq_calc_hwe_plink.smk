rule all:
    input:
        expand(
            f"reSEQ/PopGen/HWE/{{model}}/{config['species']}_{config['assembly_version']}.hardy",
            model=config["models"],
        ),


rule reseq_calc_hwe_plink:
    input:
        vcf_file=f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.vcf.gz",
        ind_file=f"Models/{config['species']}.{{model}}.ind.txt",
    output:
        f"reSEQ/PopGen/HWE/{{model}}/{config['species']}_{config['assembly_version']}.hardy",
        f"reSEQ/PopGen/HWE/{{model}}/{config['species']}_{config['assembly_version']}.log",
    log:
        f"logs/plink/{config['species']}_{config['assembly_version']}.calc_hwe_{{model}}.log",
    params:
        hwe_prefix=subpath(output[0], strip_suffix=".hardy"),
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/plink2:20240418"
    shell:
        "plink2 --vcf {input.vcf_file} --keep {input.ind_file} --hardy midp log10 --allow-extra-chr --out {params.hwe_prefix} &> {log}"
