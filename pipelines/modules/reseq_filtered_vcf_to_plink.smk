rule all:
    input:
        f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bed",
        f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bim",
        f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.fam",


rule convert_flitered_vcf_to_plink:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
    output:
        bed_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bed",
        bim_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bim",
        fam_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.fam",
    log:
        f"logs/plink/{config['species']}_{config['assembly_version']}.vcf_to_plink.log",
    params:
        out_prefix=subpath(output.bed_file, strip_suffix=".bed"),
        not_chr=(
            f"--not-chr {','.join(config['exclude_chr'])}"
            if "exclude_chr" in config and config["exclude_chr"]
            else ""
        ),
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/plink2:20240418"
    shell:
        "plink2 --vcf {input} --make-bed --out {params.out_prefix} --allow-extra-chr --set-missing-var-ids @:# --threads {threads} {params.not_chr} &> {log}"
