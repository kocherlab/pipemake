rule all:
    input:
        f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.nsl.manhattan.png",
        f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.abs_nsl.manhattan.png",


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
        mem_mb=2000,
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
        mem_mb=24000,
    threads: 12
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


rule reseq_prep_nsl_vcf_bcftools:
    input:
        "reSEQ/VCF/Phased/SplitByChrom/{chrom}.vcf.gz",
    output:
        temp("reSEQ/VCF/Phased/nSL/{chrom}.nsl.vcf.gz"),
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "bcftools view -i 'F_MISSING=0.0' {input} | bcftools annotate --set-id '%CHROM\_%POS' -O z -o {output}"


rule reseq_nsl_selscan:
    input:
        "reSEQ/VCF/Phased/nSL/{chrom}.nsl.vcf.gz",
    output:
        temp("reSEQ/PopGen/nSL/{chrom}.nsl.out"),
        temp("reSEQ/PopGen/nSL/{chrom}.nsl.log"),
    params:
        out_prefix="reSEQ/PopGen/nSL/{chrom}",
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
        "reSEQ/PopGen/nSL/{chrom}.nsl.out",
    output:
        norm_file=temp("reSEQ/PopGen/nSL/{{chrom}}.nsl.out.{config['bins']}bins.norm"),
        window_file=temp(
            "reSEQ/PopGen/nSL/{{chrom}}.nsl.out.{config['bins']}bins.norm.{str(config['window_size'])[:-3]}kb.windows"
        ),
        log_file=temp("reSEQ/PopGen/nSL/{{chrom}}.nsl.out.{config['bins']}bins.log"),
    params:
        out_prefix="reSEQ/PopGen/nSL/{chrom}.nsl.out",
        bins=config["bins"],
        window_size=config["window_size"],
    singularity:
        "docker://aewebb/selscan:v2.0.3"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        norm --nsl --bp-win --winsize {params.window_size} --files {input} --bins {params.bins} 2> {output.log_file}
        sed -i $'s/^/{wildcards.chrom}\t/' {output.window_file}
        """


def aggregate_nsl_reseq(wildcards):
    checkpoint_output = checkpoints.reseq_split_unphased_bcftools.get(
        **wildcards
    ).output[0]
    chrom_wildcards = glob_wildcards(
        os.path.join(
            checkpoint_output,
            "{chrom}.vcf.gz",
        )
    ).chrom
    return {
        "scan_nsl": expand("reSEQ/PopGen/nSL/{chrom}.nsl.out", chrom=chrom_wildcards),
        "scan_log": expand("reSEQ/PopGen/nSL/{chrom}.nsl.log", chrom=chrom_wildcards),
        "norm_nsl": expand(
            f"reSEQ/PopGen/nSL/{{chrom}}.nsl.out.{config['bins']}bins.norm",
            chrom=chrom_wildcards,
        ),
        "norm_windows": expand(
            f"reSEQ/PopGen/nSL/{{chrom}}.nsl.out.{config['bins']}bins.norm.{str(config['window_size'])[:-3]}kb.windows",
            chrom=chrom_wildcards,
        ),
        "norm_log": expand(
            f"reSEQ/PopGen/nSL/{{chrom}}.nsl.out.{config['bins']}bins.log",
            chrom=chrom_wildcards,
        ),
    }


rule reseq_cat_nsl_bash:
    input:
        unpack(aggregate_nsl_reseq),
    output:
        scan_nsl=f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.nsl.out",
        scan_log=f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.nsl.out.log",
        norm_nsl=f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.nsl.norm",
        norm_windows=f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.nsl.norm.windows",
        norm_log=f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.nsl.norm.log",
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        cat {input.scan_nsl} > {output.scan_nsl}
        cat {input.norm_nsl} > {output.norm_nsl}
        cat {input.norm_windows} > {output.norm_windows}
        cat {input.scan_log} > {output.scan_log}
        cat {input.norm_log} > {output.norm_log}
        """


rule plot_norm_nsl_pipemake:
    input:
        f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.nsl.norm",
    output:
        f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.nsl.manhattan.png",
        f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.abs_nsl.manhattan.png",
    params:
        out_prefix=f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}",
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        manhattan-plot --input-file {input} --chrom-col-int 0 --pos-col-int 0 --stat-col-int 6 --plot-stat-text "Normalized nSL" --chrom-pos-sep '_' --out-prefix {params.out_prefix}.nsl
        manhattan-plot --input-file {input} --chrom-col-int 0 --pos-col-int 0 --stat-col-int 6 --plot-stat-text "Normalized nSL" --chrom-pos-sep '_' --plot-abs --out-prefix {params.out_prefix}.abs_nsl
        """
