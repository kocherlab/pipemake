rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            f"{config['species']}_{config['assembly_version']}.xpnsl.manhattan.png",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            f"{config['species']}_{config['assembly_version']}.abs_xpnsl.manhattan.png",
        ),


rule reseq_prep_xpnsl_vcf_bcftools:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_phased_vcf_dir"],
            "SplitByChrom",
            "{chrom}.vcf.gz",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_phased_vcf_dir"],
            "SplitByChrom",
            "{chrom}.xpnsl.vcf.gz",
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


rule reseq_create_pop_xpnsl_vcf_bcftools:
    input:
        vcf=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_phased_vcf_dir"],
            "SplitByChrom",
            "{chrom}.xpnsl.vcf.gz",
        ),
        ref_inds=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["models_dir"],
            config["model_name"],
            f"{config['pop_name']}.pop",
        ),
    output:
        ref_vcf=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_phased_vcf_dir"],
            "SplitByChrom",
            "{chrom}.xpnsl.ref.vcf.gz",
        ),
        query_vcf=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_phased_vcf_dir"],
            "SplitByChrom",
            "{chrom}.xpnsl.query.vcf.gz",
        ),
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        bcftools index {input}
        bcftools view --samples-file {input.ref_inds} -O z -o {output.ref_vcf} {input.vcf}
        bcftools view --samples-file ^{input.ref_inds} -O z -o {output.query_vcf} {input.vcf}
        """


rule reseq_xpnsl_selscan:
    input:
        ref_vcf=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_phased_vcf_dir"],
            "SplitByChrom",
            "{chrom}.xpnsl.ref.vcf.gz",
        ),
        query_vcf=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_phased_vcf_dir"],
            "SplitByChrom",
            "{chrom}.xpnsl.query.vcf.gz",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "XPnSL",
                "{chrom}.xpnsl.out",
            )
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "XPnSL",
                "{chrom}.xpnsl.log",
            )
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            "{chrom}",
        ),
        maf=config["maf"],
    singularity:
        "docker://aewebb/selscan:v2.0.3"
    resources:
        mem_mb=24000,
    threads: 12
    shell:
        "selscan --xpnsl --vcf-ref {input.ref_vcf} --vcf {input.query_vcf} --maf {params.maf} --threads {threads} --out {params.out_prefix}"


rule reseq_normalize_xpnsl_norm:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            "{chrom}.xpnsl.out",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "XPnSL",
                f"{{chrom}}.xpnsl.out.norm",
            )
        ),
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "XPnSL",
                f"{{chrom}}.xpnsl.norm.log",
            )
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            "{chrom}.xpnsl",
        ),
        bins=config["bins"],
    singularity:
        "docker://aewebb/selscan:v2.0.3"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "norm --xpnsl --files {input} --bins {params.bins} 2> {params.out_prefix}.norm.log"


def aggregate_xpnsl_reseq(wildcards):
    checkpoint_output = checkpoints.reseq_split_unphased_bcftools.get(
        **wildcards
    ).output[0]
    return {
        "scan_xpnsl": expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "XPnSL",
                "{chrom}.xpnsl.out",
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
                "XPnSL",
                "{chrom}.xpnsl.log",
            ),
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.vcf.gz",
                )
            ).chrom,
        ),
        "norm_xpnsl": expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_popgen_dir"],
                "XPnSL",
                f"{{chrom}}.xpnsl.out.norm",
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
                "XPnSL",
                f"{{chrom}}.xpnsl.norm.log",
            ),
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.vcf.gz",
                )
            ).chrom,
        ),
    }


rule reseq_cat_xpnsl_bash:
    input:
        unpack(aggregate_xpnsl_reseq),
    output:
        scan_xpnsl=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            f"{config['species']}_{config['assembly_version']}.xpnsl.out",
        ),
        scan_log=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            f"{config['species']}_{config['assembly_version']}.xpnsl.out.log",
        ),
        norm_xpnsl=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            f"{config['species']}_{config['assembly_version']}.xpnsl.norm",
        ),
        norm_log=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            f"{config['species']}_{config['assembly_version']}.xpnsl.norm.log",
        ),
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        awk 'FNR>1 || NR==1' {input.scan_xpnsl} > {output.scan_xpnsl}
        awk 'FNR>1 || NR==1' {input.norm_xpnsl} > {output.norm_xpnsl}
        cat {input.scan_log} > {output.scan_log}
        cat {input.norm_log} > {output.norm_log}
        """


rule plot_norm_xpnsl_pipemake:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            f"{config['species']}_{config['assembly_version']}.xpnsl.norm",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            f"{config['species']}_{config['assembly_version']}.xpnsl.manhattan.png",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            f"{config['species']}_{config['assembly_version']}.abs_xpnsl.manhattan.png",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_popgen_dir"],
            "XPnSL",
            f"{config['species']}_{config['assembly_version']}",
        ),
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        manhattan-plot --input-file {input} --chrom-col id --pos-col id --stat-col normxpehh --plot-stat-text "Noramlized XPnSL" --chrom-pos-sep '_' --out-prefix {params.out_prefix}.xpnsl
        manhattan-plot --input-file {input} --chrom-col id --pos-col id --stat-col normxpehh --plot-stat-text "Noramlized XPnSL" --chrom-pos-sep '_' --plot-abs --out-prefix {params.out_prefix}.abs_xpnsl
        """
