rule all:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",


rule filter_basic_vcf_bcftools:
    input:
        vcf_file=f"reSEQ/VCF/Unfiltered/{config['species']}_{config['assembly_version']}.vcf.gz",
        ind_file=f"Models/{config['species']}.{config['model_name']}.ind.txt",
    output:
        vcf_file=temp(
            f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.woMD.vcf.gz"
        ),
        idx_file=temp(
            f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.woMD.vcf.gz.csi"
        ),
    log:
        f"logs/bcftools/{config['species']}_{config['assembly_version']}.filter.log",
    params:
        min_alleles=config["min_alleles"],
        max_alleles=config["max_alleles"],
        maf=config["maf_cutoff"],
        qual=config["qual_cutoff"],
    resources:
        mem_mb=16000,
    threads: 4
    singularity:
        "docker://aewebb/bcftools:v1.20"
    shell:
        """
        bcftools view --samples-file {input.ind_file} {input.vcf_file} 2> {log} | bcftools view --min-alleles {params.min_alleles} --max-alleles {params.max_alleles} --types snps --include 'MAF>={params.maf} && QUAL>={params.qual}' --output-type z --output-file {output.vcf_file} --threads {threads} 2>> {log}
        bcftools index {output.vcf_file} 2>> {log}
        """


checkpoint pop_ind_file:
    input:
        f"Models/{config['species']}.model",
    output:
        temp(directory("Models/BCFtools")),
    log:
        f"logs/model-pop-files/{config['species']}.log",
    params:
        model_name=config["model_name"],
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.4.3"
    shell:
        "model-pop-files --model-file {input} --model-name {params.model_name} --out-dir {output} &> {log}"


rule pop_vcf_bcftools:
    input:
        vcf_file=f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.woMD.vcf.gz",
        pop_file="Models/BCFtools/{model_pop}.pop",
    output:
        vcf_file=temp(
            f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.woMD.{{model_pop}}.vcf.gz"
        ),
        idx_file=temp(
            f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.woMD.{{model_pop}}.vcf.gz.csi"
        ),
    log:
        f"logs/bcftools/{{model_pop}}.filter_pop.log",
    params:
        missing_cutoff=config["missing_cutoff"],
    resources:
        mem_mb=8000,
    threads: 1
    singularity:
        "docker://aewebb/bcftools:v1.20"
    shell:
        """
        bcftools view --samples-file {input.pop_file} {input.vcf_file} 2> {log} | bcftools view -i 'F_MISSING<={params.missing_cutoff}' --output-type z --output-file {output.vcf_file} 2>> {log}
        bcftools index {output.vcf_file} 2>> {log}
        """


def aggregate_pop_reseq(wildcards):
    checkpoint_output = checkpoints.pop_ind_file.get(**wildcards).output[0]
    return {
        "vcf_file": expand(
            f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.woMD.{{model_pop}}.vcf.gz",
            model_pop=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{model_pop}.pop",
                )
            ).model_pop,
        ),
        "idx_file": expand(
            f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.woMD.{{model_pop}}.vcf.gz.csi",
            model_pop=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{model_pop}.pop",
                )
            ).model_pop,
        ),
    }


rule isec_pop_vcfs_bcftools:
    input:
        unpack(aggregate_pop_reseq),
    output:
        temp(
            f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.missing_data.sites"
        ),
    log:
        f"logs/bcftools/{config['species']}_{config['assembly_version']}.isec.log",
    params:
        pop_count=lambda wildcards, input: len(input) / 2,
    resources:
        mem_mb=16000,
    threads: 4
    singularity:
        "docker://aewebb/bcftools:v1.20"
    shell:
        "bcftools isec {input.vcf_file} -n={params.pop_count} 2> {log} | cut -f1,2 > {output}"


rule filter_pops_missing_data_vcf_bcftools:
    input:
        vcf_file=f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.woMD.vcf.gz",
        idx_file=f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.woMD.vcf.gz.csi",
        sites_file=f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.missing_data.sites",
    output:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
    log:
        f"logs/bcftools/{config['species']}_{config['assembly_version']}.filter_missing_data.log",
    resources:
        mem_mb=8000,
    threads: 4
    singularity:
        "docker://aewebb/bcftools:v1.20"
    shell:
        "bcftools view -R {input.sites_file} --output-type z --output-file {output} --threads {threads} {input.vcf_file} 2> {log}"
