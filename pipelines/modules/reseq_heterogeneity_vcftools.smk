rule all:
    input:
        f"reSEQ/PopGen/vcftools/{config['species']}_{config['assembly_version']}.sample.het",


rule het_region_vcftools:
    input:
        f"reSEQ/VCF/{config['species']}_{config['assembly_version']}.vcf.gz",
    output:
        f"reSEQ/PopGen/vcftools/{config['species']}_{config['assembly_version']}.site.het",
    log:
        f"logs/vcftools/{config['species']}_{config['assembly_version']}.site.het.log",
    params:
        out_prefix=subpath(output[0], strip_suffix=".het"),
        chrom=f"--chr {config['chrom']}" if "chrom" in config else "",
        start=f"--from-bp {config['start']}" if "start" in config else "",
        end=f"--to-bp {config['end']}" if "end" in config else "",
    container:
        "docker://aewebb/vcftools:v0.1.17"
    shell:
        "vcftools --gzvcf {input} {params.chrom} {params.start} {params.end} --het --out {params.out_prefix} &> {log}"


rule het_per_indv:
    input:
        f"reSEQ/PopGen/vcftools/{config['species']}_{config['assembly_version']}.site.het",
    output:
        f"reSEQ/PopGen/vcftools/{config['species']}_{config['assembly_version']}.sample.het2",
    shell:
        """awk 'NR==1{{print "INDV\\tN_SITES\\tN_HET\\tHET_FRAC\\tF"; next}} {{het=$4-$2; printf "%s\\t%d\\t%d\\t%.6f\\t%s\\n",$1,$4,het,het/$4,$5}}' {input} > {output}"""
