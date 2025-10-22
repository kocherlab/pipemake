rule all:
    input:
        f"reSEQ/PopGen/XPnSL/{config['species']}_{config['assembly_version']}.xpnsl.manhattan.png",
        f"reSEQ/PopGen/XPnSL/{config['species']}_{config['assembly_version']}.abs_xpnsl.manhattan.png",


checkpoint reseq_split_unphased_bcftools:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.vcf.gz",
    output:
        temp(directory("reSEQ/VCF/Filtered/SplitByChrom")),
    params:
        out_dir="reSEQ/VCF/Filtered/SplitByChrom",
        out_prefix="reSEQ/VCF/Filtered/SplitByChrom/",
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        """
        mkdir {params.out_dir}
        bcftools index -f {input}
        bcftools index -s {input} | cut -f 1 | while read chrom; do bcftools view --regions $chrom -O z -o {params.out_prefix}${{chrom}}.vcf.gz {input}; done
        """


rule reseq_index_unphased_bcftools:
    input:
        "reSEQ/VCF/Filtered/SplitByChrom/{chrom}.vcf.gz",
    output:
        "reSEQ/VCF/Filtered/SplitByChrom/{chrom}.vcf.gz.csi",
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "bcftools index -f {input}"


rule reseq_phase_chroms_shapeit4:
    input:
        vcf="reSEQ/VCF/Filtered/SplitByChrom/{chrom}.vcf.gz",
        index="reSEQ/VCF/Filtered/SplitByChrom/{chrom}.vcf.gz.csi",
    output:
        temp("reSEQ/VCF/Phased/SplitByChrom/{chrom}.vcf.gz"),
    singularity:
        "docker://aewebb/shapeit4:v4.2.2"
    resources:
        mem_mb=24000,
    threads: 12
    shell:
        "shapeit4 --input {input.vcf} --region {wildcards.chrom} --output {output} --thread {threads}"


rule reseq_create_pop_xpnsl_vcf_bcftools:
    input:
        vcf="reSEQ/VCF/Phased/SplitByChrom/{chrom}.vcf.gz",
        ref_inds=f"Models/{config['model_name']}/{config['pop_name']}.pop",
    output:
        ref_vcf="reSEQ/VCF/Phased/SplitByChrom/{chrom}.ref.vcf.gz",
        query_vcf="reSEQ/VCF/Phased/SplitByChrom/{chrom}.query.vcf.gz",
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
        ref_vcf="reSEQ/VCF/Phased/SplitByChrom/{chrom}.ref.vcf.gz",
        query_vcf="reSEQ/VCF/Phased/SplitByChrom/{chrom}.query.vcf.gz",
    output:
        temp("reSEQ/PopGen/XPnSL/{chrom}.xpnsl.out"),
        temp("reSEQ/PopGen/XPnSL/{chrom}.xpnsl.log"),
    params:
        out_prefix="reSEQ/PopGen/XPnSL/{chrom}",
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
        "reSEQ/PopGen/XPnSL/{chrom}.xpnsl.out",
    output:
        temp("reSEQ/PopGen/XPnSL/{{chrom}}.xpnsl.out.norm"),
        temp("reSEQ/PopGen/XPnSL/{{chrom}}.xpnsl.norm.log"),
    params:
        out_prefix="reSEQ/PopGen/XPnSL/{chrom}.xpnsl",
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
            "reSEQ/PopGen/XPnSL/{chrom}.xpnsl.out",
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.vcf.gz",
                )
            ).chrom,
        ),
        "scan_log": expand(
            "reSEQ/PopGen/XPnSL/{chrom}.xpnsl.log",
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.vcf.gz",
                )
            ).chrom,
        ),
        "norm_xpnsl": expand(
            "reSEQ/PopGen/XPnSL/{chrom}.xpnsl.out.norm",
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.vcf.gz",
                )
            ).chrom,
        ),
        "norm_log": expand(
            "reSEQ/PopGen/XPnSL/{chrom}.xpnsl.norm.log",
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
        scan_xpnsl=f"reSEQ/PopGen/XPnSL/{config['species']}_{config['assembly_version']}.xpnsl.out",
        scan_log=f"reSEQ/PopGen/XPnSL/{config['species']}_{config['assembly_version']}.xpnsl.out.log",
        norm_xpnsl=f"reSEQ/PopGen/XPnSL/{config['species']}_{config['assembly_version']}.xpnsl.norm",
        norm_log=f"reSEQ/PopGen/XPnSL/{config['species']}_{config['assembly_version']}.xpnsl.norm.log",
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
        f"reSEQ/PopGen/XPnSL/{config['species']}_{config['assembly_version']}.xpnsl.norm",
    output:
        f"reSEQ/PopGen/XPnSL/{config['species']}_{config['assembly_version']}.xpnsl.manhattan.png",
        f"reSEQ/PopGen/XPnSL/{config['species']}_{config['assembly_version']}.abs_xpnsl.manhattan.png",
    params:
        out_prefix=f"reSEQ/PopGen/XPnSL/{config['species']}_{config['assembly_version']}",
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        """
        manhattan-plot --input-file {input} --chrom-col id --pos-col id --stat-col normxpehh --plot-stat-text "Noramlized XPnSL" --chrom-pos-sep '_' --out-prefix {params.out_prefix}.xpnsl
        manhattan-plot --input-file {input} --chrom-col id --pos-col id --stat-col normxpehh --plot-stat-text "Noramlized XPnSL" --chrom-pos-sep '_' --plot-abs --out-prefix {params.out_prefix}.abs_xpnsl
        """
