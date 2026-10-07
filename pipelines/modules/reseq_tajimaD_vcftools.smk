rule all:
    input:
        f"reSEQ/PopGen/vcftools/{config['species']}_{config['assembly_version']}.Tajima.D",


rule tajimaD_vcftools:
    input:
        f"reSEQ/VCF/{config['species']}_{config['assembly_version']}.vcf.gz",
    output:
        f"reSEQ/PopGen/vcftools/{config['species']}_{config['assembly_version']}.Tajima.D",
    log:
        f"logs/vcftools/{config['species']}_{config['assembly_version']}.Tajima.D.log",
    params:
        out_prefix=subpath(output[0], strip_suffix=".Tajima.D"),
        chrom=f"--chr {config['chrom']}" if "chrom" in config else "",
        start=f"--from-bp {config['start']}" if "start" in config else "",
        end=f"--to-bp {config['end']}" if "end" in config else "",
        window=config["tajima_window"],
    container:
        "docker://aewebb/vcftools:v0.1.17"
    shell:
        "vcftools --gzvcf {input} {params.chrom} {params.start} {params.end} --TajimaD {params.window} --out {params.out_prefix} &> {log}"
