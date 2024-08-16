rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_filtered_vcf_dir"],
            f"{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
        ),


rule filter_basic_vcf_bcftools:
    input:
        vcf_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_unfiltered_vcf_dir"],
            f"{config['species']}_{config['assembly_version']}.vcf.gz",
        ),
        ind_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.{config['model_name']}.ind.txt",
        ),
    output:
        vcf_file=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_filtered_vcf_dir"],
                f"{config['species']}_{config['assembly_version']}.filtered.woMD.vcf.gz",
            )
        ),
        idx_file=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_filtered_vcf_dir"],
                f"{config['species']}_{config['assembly_version']}.filtered.woMD.vcf.gz.csi",
            )
        ),
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
        bcftools view --samples-file {input.ind_file} {input.vcf_file} | bcftools view --min-alleles {params.min_alleles} --max-alleles {params.max_alleles} --types snps --include 'MAF>={params.maf} && QUAL>={params.qual}' --output-type z --output-file {output.vcf_file} --threads {threads}
        bcftools index {output.vcf_file}
        """


checkpoint pop_ind_file:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            f"{config['species']}.model",
        ),
    output:
        temp(
            directory(
                os.path.join(
                    config["paths"]["workflow_prefix"],
                    config["paths"]["models_dir"],
                    "BCFtools",
                )
            )
        ),
    params:
        model_name=config["model_name"],
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v0.1.27"
    shell:
        "model-pop-files --model-file {input} --model-name {params.model_name} --out-dir {output}"


rule pop_vcf_bcftools:
    input:
        vcf_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_filtered_vcf_dir"],
            f"{config['species']}_{config['assembly_version']}.filtered.woMD.vcf.gz",
        ),
        pop_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            "BCFtools",
            "{model_pop}.pop",
        ),
    output:
        vcf_file=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_filtered_vcf_dir"],
                f"{config['species']}_{config['assembly_version']}.filtered.woMD.{{model_pop}}.vcf.gz",
            )
        ),
        idx_file=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_filtered_vcf_dir"],
                f"{config['species']}_{config['assembly_version']}.filtered.woMD.{{model_pop}}.vcf.gz.csi",
            )
        ),
    params:
        missing_cutoff=config["missing_cutoff"],
    resources:
        mem_mb=8000,
    threads: 1
    singularity:
        "docker://aewebb/bcftools:v1.20"
    shell:
        """
        bcftools view --samples-file {input.pop_file} {input.vcf_file} | bcftools view -i 'F_MISSING<={params.missing_cutoff}' --output-type z --output-file {output.vcf_file}
        bcftools index {output.vcf_file}
        """


def aggregate_pop_reseq(wildcards):
    checkpoint_output = checkpoints.pop_ind_file.get(**wildcards).output[0]
    return {
        "vcf_file": expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_filtered_vcf_dir"],
                f"{config['species']}_{config['assembly_version']}.filtered.woMD.{{model_pop}}.vcf.gz",
            ),
            model_pop=glob_wildcards(
                os.path.join(
                    config["paths"]["workflow_prefix"],
                    checkpoint_output,
                    "{model_pop}.pop",
                )
            ).model_pop,
        ),
        "idx_file": expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_filtered_vcf_dir"],
                f"{config['species']}_{config['assembly_version']}.filtered.woMD.{{model_pop}}.vcf.gz.csi",
            ),
            model_pop=glob_wildcards(
                os.path.join(
                    config["paths"]["workflow_prefix"],
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
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_filtered_vcf_dir"],
                f"{config['species']}_{config['assembly_version']}.filtered.missing_data.sites",
            )
        ),
    params:
        pop_count=lambda wildcards, input: len(input) / 2,
    resources:
        mem_mb=16000,
    threads: 4
    singularity:
        "docker://aewebb/bcftools:v1.20"
    shell:
        "bcftools isec {input.vcf_file} -n={params.pop_count} | cut -f1,2 > {output}"


rule filter_pops_missing_data_vcf_bcftools:
    input:
        vcf_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_filtered_vcf_dir"],
            f"{config['species']}_{config['assembly_version']}.filtered.woMD.vcf.gz",
        ),
        idx_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_filtered_vcf_dir"],
            f"{config['species']}_{config['assembly_version']}.filtered.woMD.vcf.gz.csi",
        ),
        sites_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_filtered_vcf_dir"],
            f"{config['species']}_{config['assembly_version']}.filtered.missing_data.sites",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_filtered_vcf_dir"],
            f"{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
        ),
    resources:
        mem_mb=8000,
    threads: 4
    singularity:
        "/Genomics/kocherlab/lab/Pipelines/images/kocherPOP.sif"
    shell:
        "bcftools view -R {input.sites_file} --output-type z --output-file {output} --threads {threads} {input.vcf_file}"
