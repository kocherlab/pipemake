rule all:
    input:
        expand("ATAC_seq/MACS3/{sample}_peaks.narrowPeak", sample=config["samples"]),


ruleorder: MACS3_peaks_pair_end > MACS3_peaks_single_end


rule MACS3_peaks_pair_end:
    input:
        bam="ATAC_seq/BAM/Deduplicated/{sample}.deduplicated.bam",
        r1_reads="ATAC_seq/FASTQ/{sample}_R1.fastq.gz",
        r2_reads="ATAC_seq/FASTQ/{sample}_R2.fastq.gz",
    output:
        "ATAC_seq/MACS3/{sample}_peaks.narrowPeak",
    log:
        "logs/macs3/{sample}.peaks.log",
    params:
        out_prefix=subpath(output[0], strip_suffix="_peaks.narrowPeak"),
        gs=config["genome_size"],
    singularity:
        "docker://aewebb/macs3:v3.0.1"
    resources:
        mem_mb=64000,
    threads: 4
    shell:
        "macs3 callpeak -t {input.bam} -f BAMPE -g {params.gs} -n {params.out_prefix} --keep-dup all &> {log}"


rule MACS3_peaks_single_end:
    input:
        bam="ATAC_seq/BAM/Deduplicated/{sample}.deduplicated.bam",
        r1_reads="ATAC_seq/FASTQ/{sample}_R1.fastq.gz",
    output:
        "ATAC_seq/MACS3/{sample}_peaks.narrowPeak",
    log:
        "logs/macs3/{sample}.peaks.log",
    params:
        out_prefix=subpath(output[0], strip_suffix="_peaks.narrowPeak"),
        gs=config["genome_size"],
    singularity:
        "docker://aewebb/macs3:v3.0.1"
    resources:
        mem_mb=64000,
    threads: 4
    shell:
        "macs3 callpeak -t {input.bam} -f BAM -g {params.gs} -n {params.out_prefix} --keep-dup all &> {log}"
