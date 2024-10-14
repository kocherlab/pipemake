rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            f"{config['species']}_{config['assembly_version']}.nsl.manhattan.png",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            f"{config['species']}_{config['assembly_version']}.abs_nsl.manhattan.png",
        ),


rule reseq_prep_nsl_vcf_bcftools:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_phased_vcf_dir"],
            "SplitByChrom",
            "{chrom}.vcf.gz",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_phased_vcf_dir"],
                "SplitByChrom",
                "{chrom}.nsl.vcf.gz",
            )
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_phased_vcf_dir"],
            "SplitByChrom",
            "{chrom}",
        ),
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "bcftools view -i 'F_MISSING=0.0' {input} | bcftools annotate --set-id '%CHROM\_%POS' -O z -o {output}"


rule reseq_nsl_selscan:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_phased_vcf_dir"],
            "SplitByChrom",
            "{chrom}.nsl.vcf.gz",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "nSL",
                "{chrom}.nsl.out",
            )
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "nSL",
                "{chrom}.nsl.log",
            )
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            "{chrom}",
        ),
        maf=config["maf"],
    singularity:
        "docker://aewebb/selscan:v2.0.3"
    resources:
        mem_mb=24000,
    threads: 12
    shell:
        "selscan --nsl --vcf {input} --maf {params.maf} --threads {threads} --out {params.out_prefix}"


rule reseq_normalize_nsl_norm:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            "{chrom}.nsl.out",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "nSL",
                f"{{chrom}}.nsl.out.{config['bins']}bins.norm",
            )
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "nSL",
                f"{{chrom}}.nsl.out.{config['bins']}bins.log",
            )
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            "{chrom}.nsl.out",
        ),
        bins=config["bins"],
    singularity:
        "docker://aewebb/selscan:v2.0.3"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "norm --nsl --files {input} --bins {params.bins} 2> {params.out_prefix}.{params.bins}bins.log"


def aggregate_nsl_reseq(wildcards):
    checkpoint_output = checkpoints.reseq_split_unphased_bcftools.get(
        **wildcards
    ).output[0]
    return {
        "scan_nsl": expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "nSL",
                "{chrom}.nsl.out",
            ),
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.vcf.gz",
                )
            ).chrom,
        ),
        "scan_log": expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "nSL",
                "{chrom}.nsl.log",
            ),
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.vcf.gz",
                )
            ).chrom,
        ),
        "norm_nsl": expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "nSL",
                f"{{chrom}}.nsl.out.{config['bins']}bins.norm",
            ),
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.vcf.gz",
                )
            ).chrom,
        ),
        "norm_log": expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "nSL",
                f"{{chrom}}.nsl.out.{config['bins']}bins.log",
            ),
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.vcf.gz",
                )
            ).chrom,
        ),
    }


rule reseq_cat_nsl_bash:
    input:
        unpack(aggregate_nsl_reseq),
    output:
        scan_nsl=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            f"{config['species']}_{config['assembly_version']}.nsl.out",
        ),
        scan_log=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            f"{config['species']}_{config['assembly_version']}.nsl.out.log",
        ),
        norm_nsl=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            f"{config['species']}_{config['assembly_version']}.nsl.norm",
        ),
        norm_log=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            f"{config['species']}_{config['assembly_version']}.nsl.norm.log",
        ),
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        cat {input.scan_nsl} > {output.scan_nsl}
        cat {input.norm_nsl} > {output.norm_nsl}
        cat {input.scan_log} > {output.scan_log}
        cat {input.norm_log} > {output.norm_log}
        """


rule plot_norm_nsl_pipemake:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            f"{config['species']}_{config['assembly_version']}.nsl.norm",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            f"{config['species']}_{config['assembly_version']}.nsl.manhattan.png",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            f"{config['species']}_{config['assembly_version']}.abs_nsl.manhattan.png",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "nSL",
            f"{config['species']}_{config['assembly_version']}",
        ),
    singularity:
        "docker://aewebb/pipemake_utils:v0.1.27"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        manhattan-plot --input-file {input} --chrom-col-int 0 --pos-col-int 0 --stat-col-int 6 --plot-stat-text "Noramlized nSL" --chrom-pos-sep '_' --out-prefix {params.out_prefix}.nsl
        manhattan-plot --input-file {input} --chrom-col-int 0 --pos-col-int 0 --stat-col-int 6 --plot-stat-text "Noramlized nSL" --chrom-pos-sep '_' --plot-abs --out-prefix {params.out_prefix}.abs_nsl
        """
