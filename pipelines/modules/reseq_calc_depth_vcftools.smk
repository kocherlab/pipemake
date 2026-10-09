rule all:
    input:
        f"reSEQ/PopGen/vcftools/{config['species']}_{config['assembly_version']}.idepth",


rule depth_per_indv_vcftools:
    input:
        f"reSEQ/VCF/{config['species']}_{config['assembly_version']}.vcf.gz",
    output:
        f"reSEQ/PopGen/vcftools/{config['species']}_{config['assembly_version']}.idepth",
    log:
        f"logs/vcftools/{config['species']}_{config['assembly_version']}.idepth.log",
    params:
        out_prefix=subpath(output[0], strip_suffix=".idepth"),
        chrom=f"--chr {config['chrom']}" if "chrom" in config else "",
        start=f"--from-bp {config['start']}" if "start" in config else "",
        end=f"--to-bp {config['end']}" if "end" in config else "",
    container:
        "docker://aewebb/vcftools:v0.1.17"
    shell:
        "vcftools --gzvcf {input} {params.chrom} {params.start} {params.end} --depth --out {params.out_prefix} &> {log}"
