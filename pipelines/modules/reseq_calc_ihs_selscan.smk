rule all:
    input:
        f"reSEQ/PopGen/ihs/{config['species']}_{config['assembly_version']}.ihs.out",
        f"reSEQ/PopGen/ihs/{config['species']}_{config['assembly_version']}.ihs.out.log",
        f"reSEQ/PopGen/ihs/{config['species']}_{config['assembly_version']}.ihs.norm",
        f"reSEQ/PopGen/ihs/{config['species']}_{config['assembly_version']}.ihs.norm.log",


rule reseq_phased_map_plink:
    input:
        "reSEQ/VCF/Phased/SplitByChrom/{sweep_chrom}.vcf.gz",
    output:
        temp("reSEQ/VCF/Phased/SplitByChrom/{sweep_chrom}.map"),
        temp("reSEQ/VCF/Phased/SplitByChrom/{sweep_chrom}.ped"),
        temp("reSEQ/VCF/Phased/SplitByChrom/{sweep_chrom}.log"),
    params:
        out_prefix="reSEQ/VCF/Phased/SplitByChrom/{sweep_chrom}",
    singularity:
        "docker://aewebb/plink2:20240418"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "plink2 --vcf {input} --export ped --out {params.out_prefix} --set-all-var-ids @:# --allow-extra-chr"


rule reseq_phased_ids_bcftools:
    input:
        "reSEQ/VCF/Phased/SplitByChrom/{sweep_chrom}.vcf.gz",
    output:
        temp("reSEQ/VCF/Phased/SplitByChrom/{sweep_chrom}.id.vcf.gz"),
    params:
        out_prefix="reSEQ/VCF/Phased/SplitByChrom/{sweep_chrom}",
    singularity:
        "docker://aewebb/bcftools:v1.20"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "bcftools annotate --set-id '%CHROM\_%POS' {input} -O z -o {output}"


rule reseq_ihs_selscan:
    input:
        vcf="reSEQ/VCF/Phased/SplitByChrom/{sweep_chrom}.id.vcf.gz",
        map="reSEQ/VCF/Phased/SplitByChrom/{sweep_chrom}.map",
    output:
        temp("reSEQ/PopGen/ihs/{sweep_chrom}.ihs.out"),
        temp("reSEQ/PopGen/ihs/{sweep_chrom}.ihs.log"),
    params:
        out_prefix="reSEQ/PopGen/ihs/{sweep_chrom}",
        maf=config["maf"],
    singularity:
        "docker://aewebb/selscan:v2.0.3"
    resources:
        mem_mb=24000,
    threads: 12
    shell:
        "selscan --ihs --ihs-detail --vcf {input.vcf} --map {input.map} --pmap --maf {params.maf} --threads {threads} --out {params.out_prefix}"


rule reseq_ihs_normalize_norm:
    input:
        "reSEQ/PopGen/ihs/{sweep_chrom}.ihs.out",
    output:
        temp(f"reSEQ/PopGen/ihs/{{sweep_chrom}}.ihs.out.{config['bins']}bins.norm"),
        temp(f"reSEQ/PopGen/ihs/{{sweep_chrom}}.ihs.out.{config['bins']}bins.log"),
    params:
        out_prefix="reSEQ/PopGen/ihs/{sweep_chrom}.ihs.out",
        bins=config["bins"],
    singularity:
        "docker://aewebb/selscan:v2.0.3"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "norm --ihs --files {input} --bins {params.bins} 2> {params.out_prefix}.{params.bins}bins.log"


def aggregate_ihs_reseq(wildcards):
    checkpoint_output = checkpoints.reseq_split_unphased_bcftools.get(
        **wildcards
    ).output[0]
    return {
        "scan_ihs": expand(
            "reSEQ/PopGen/ihs/{sweep_chrom}.ihs.out",
            sweep_chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{sweep_chrom}.vcf.gz",
                )
            ).sweep_chrom,
        ),
        "scan_log": expand(
            "reSEQ/PopGen/ihs/{sweep_chrom}.ihs.log",
            sweep_chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{sweep_chrom}.vcf.gz",
                )
            ).sweep_chrom,
        ),
        "norm_ihs": expand(
            f"reSEQ/PopGen/hs/{{sweep_chrom}}.ihs.out.{config['bins']}bins.norm",
            sweep_chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{sweep_chrom}.vcf.gz",
                )
            ).sweep_chrom,
        ),
        "norm_log": expand(
            f"reSEQ/PopGen/ihs/{{sweep_chrom}}.ihs.out.{config['bins']}bins.log",
            sweep_chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{sweep_chrom}.vcf.gz",
                )
            ).sweep_chrom,
        ),
    }


rule reseq_cat_ihs_bash:
    input:
        unpack(aggregate_ihs_reseq),
    output:
        scan_ihs=f"reSEQ/PopGen/ihs/{config['species']}_{config['assembly_version']}.ihs.out",
        scan_log=f"reSEQ/PopGen/ihs/{config['species']}_{config['assembly_version']}.ihs.out.log",
        norm_ihs=f"reSEQ/PopGen/ihs/{config['species']}_{config['assembly_version']}.ihs.norm",
        norm_log=f"reSEQ/PopGen/ihs/{config['species']}_{config['assembly_version']}.ihs.norm.log",
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        cat {input.scan_ihs} > {output.scan_ihs}
        cat {input.norm_ihs} > {output.norm_ihs}
        cat {input.scan_log} > {output.scan_log}
        cat {input.norm_log} > {output.norm_log}
        """
