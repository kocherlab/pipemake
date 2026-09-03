rule all:
    input:
        expand("FASTQ/Filtered/{sample}.filt.fastq.gz", sample=config["samples"]),


rule hifi_reads_screen_hifiadapterfilt:
    input:
        "FASTQ/Unfiltered/{sample}.fastq.gz",
    output:
        "FASTQ/Filtered/{sample}.filt.fastq.gz",
        temp("FASTQ/Filtered/{sample}.blocklist"),
        temp("FASTQ/Filtered/{sample}.contaminant.blastout"),
    log:
        "logs/hifiadapterfilt/{sample}.log",
    params:
        input_prefix=subpath(input[0], strip_suffix=".fastq.gz"),
        output_dir=subpath(output[0], parent=True),
        min_length=config["min_length"],
        min_match=config["min_match"],
    singularity:
        "docker://aewebb/hifiadapterfilt:v3.0.0"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "hifiadapterfilt.sh -p {params.input_prefix} -l {params.min_length} -m {params.min_match} -o {params.output_dir} -t {threads} &> {log}"
