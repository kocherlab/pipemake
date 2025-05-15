rule all:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",


rule convert_flitered_vcf_to_plink:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
    output:
        f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bed",
        f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bim",
        f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.fam",
    params:
        out_prefix=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered",
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/plink2:20240418"
    shell:
        "plink2 --vcf {input} --make-bed --out {params.out_prefix} --allow-extra-chr --set-missing-var-ids @:# --threads {threads}"
