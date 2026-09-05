rule all:
    input:
        f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.nsl.manhattan.png",
        f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}.abs_nsl.manhattan.png",


checkpoint reseq_split_unphased_bcftools:
    input:
        f"reSEQ/VCF/Filtered/{config['species']}_{config['assembly_version']}.vcf.gz",
    output:
        temp(directory("reSEQ/VCF/Filtered/SplitByChrom")),
    log:
        f"logs/bcftools/{config['species']}_{config['assembly_version']}.split_by_chrom.log",
    params:
        out_prefix=lambda wildcards, output: os.path.join(output[0],),
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        """
        mkdir {output}
        bcftools index -f {input}
        bcftools index -s {input} | cut -f 1 | while read chrom; do bcftools view --regions $chrom -O z -o {params.out_prefix}${{chrom}}.vcf.gz {input} 2> {log}; done
        """


rule reseq_index_unphased_bcftools:
    input:
        "reSEQ/VCF/Filtered/SplitByChrom/{chrom}.vcf.gz",
    output:
        "reSEQ/VCF/Filtered/SplitByChrom/{chrom}.vcf.gz.csi",
    log:
        "logs/bcftools/{chrom}.index.log",
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "bcftools index -f {input} 2> {log}"


rule reseq_phase_chroms_shapeit4:
    input:
        vcf="reSEQ/VCF/Filtered/SplitByChrom/{chrom}.vcf.gz",
        index="reSEQ/VCF/Filtered/SplitByChrom/{chrom}.vcf.gz.csi",
    output:
        temp("reSEQ/VCF/Phased/SplitByChrom/{chrom}.vcf.gz"),
    log:
        "logs/shapeit4/{chrom}.phase.log",
    singularity:
        "docker://aewebb/shapeit4:v4.2.2"
    resources:
        mem_mb=24000,
    threads: 12
    shell:
        "shapeit4 --input {input.vcf} --region {wildcards.chrom} --output {output} --thread {threads} --log {log}"


rule reseq_remove_missing_data_chroms_bcftools:
    input:
        "reSEQ/VCF/Phased/SplitByChrom/{chrom}.vcf.gz",
    output:
        temp("reSEQ/VCF/nSL/SplitByChrom/{chrom}.vcf.gz"),
    log:
        "logs/bcftools/{chrom}.remove_missing.log",
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "bcftools view -i 'F_MISSING=0.0' {input} -O z -o {output} 2> {log}"


rule reseq_nsl_selscan:
    input:
        "reSEQ/VCF/nSL/SplitByChrom/{chrom}.vcf.gz",
    output:
        temp("reSEQ/PopGen/nSL/{chrom}.nsl.out"),
    log:
        "logs/selscan/{chrom}.nsl.log",
    params:
        out_prefix=subpath(output[0], strip_suffix=".nsl.out"),
        maf=config["maf"],
    singularity:
        "docker://aewebb/selscan:v2.0.3"
    resources:
        mem_mb=24000,
    threads: 12
    shell:
        """
        selscan --nsl --vcf {input} --maf {params.maf} --threads {threads} --out {params.out_prefix}
        mv {params.out_prefix}.nsl.log {log}
        """


rule reseq_normalize_nsl_norm:
    input:
        "reSEQ/PopGen/nSL/{chrom}.nsl.out",
    output:
        norm_file=temp(f"reSEQ/PopGen/nSL/{{chrom}}.nsl.out.{config['bins']}bins.norm"),
        window_file=temp(
            f"reSEQ/PopGen/nSL/{{chrom}}.nsl.out.{config['bins']}bins.norm.{str(config['window_size'])[:-3]}kb.windows"
        ),
    log:
        f"logs/selscan/{{chrom}}.nsl.out.{config['bins']}bins.log",
    params:
        bins=config["bins"],
        window_size=config["window_size"],
    singularity:
        "docker://aewebb/selscan:v2.0.3"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        norm --nsl --bp-win --winsize {params.window_size} --files {input} --bins {params.bins} 2> {log}
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
    log:
        f"logs/plot-nsl/{config['species']}_{config['assembly_version']}.log",
    params:
        out_prefix=f"reSEQ/PopGen/nSL/{config['species']}_{config['assembly_version']}",
    singularity:
        "docker://aewebb/pipemake_utils:v1.4.3"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        """
        manhattan-plot --input-file {input} --chrom-col-int 0 --pos-col-int 0 --stat-col-int 6 --plot-stat-text "Normalized nSL" --chrom-pos-sep '_' --out-prefix {params.out_prefix}.nsl &> {log}
        manhattan-plot --input-file {input} --chrom-col-int 0 --pos-col-int 0 --stat-col-int 6 --plot-stat-text "Normalized nSL" --chrom-pos-sep '_' --plot-abs --out-prefix {params.out_prefix}.abs_nsl &>> {log}
        """
