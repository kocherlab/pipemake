rule all:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.pops.log",


checkpoint pop_ind_file:
    input:
        f"Models/{config['species']}.model",
    output:
        temp(directory("Models/BCFtools")),
    params:
        model_name=config["model_name"],
    resources:
        mem_mb=2000,
    threads: 1
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    shell:
        "model-pop-files --model-file {input} --model-name {params.model_name} --out-dir {output}"


rule pop_vcf_bcftools:
    input:
        vcf_file=f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
        pop_file="Models/BCFtools/{model_pop}.pop",
    output:
        f"reSEQ/VCF/Filtered/{{model_pop}}/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
    resources:
        mem_mb=8000,
    threads: 1
    singularity:
        "docker://aewebb/bcftools:v1.20"
    shell:
        "bcftools view --samples-file {input.pop_file} --output-type z --output-file {output} {input.vcf_file}"


def aggregate_pop_reseq(wildcards):
    checkpoint_output = checkpoints.pop_ind_file.get(**wildcards).output[0]
    return expand(
        f"reSEQ/VCF/Filtered/{{model_pop}}/{config['species']}_{config['assembly_version']}.filtered.vcf.gz",
        model_pop=glob_wildcards(
            os.path.join(
                checkpoint_output,
                "{model_pop}.pop",
            )
        ).model_pop,
    )


rule log_pop_vcfs_bash:
    input:
        unpack(aggregate_pop_reseq),
    output:
        temp(
            f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.pops.log"
        ),
    params:
        pop_count=lambda wildcards, input: len(input) / 2,
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "echo {input} > {output}"
