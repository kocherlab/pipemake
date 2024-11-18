rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["atacseq_peak_dir"],
                "MACS3",
                "{sample}_peaks.narrowPeak",
            ),
            sample=config["samples"],
        ),


ruleorder: MACS3_peaks_pair_end > MACS3_peaks_single_end


rule MACS3_peaks_pair_end:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_dedup_bam_dir"],
            "{sample}.deduplicated.bam",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_peak_dir"],
            "MACS3",
            "{sample}_peaks.narrowPeak",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_peak_dir"],
            "MACS3",
            "{sample}",
        ),
        gs=config["genome_size"],
    singularity:
        "docker://aewebb/macs3:v3.0.1"
    resources:
        mem_mb=4000,
    threads: 4
    shell:
        "macs3 callpeak -t {input} -f BAMPE -g {params.gs} -n {params.out_prefix} --keep-dup all"


rule MACS3_peaks_single_end:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_dedup_bam_dir"],
            "{sample}.deduplicated.bam",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_peak_dir"],
            "MACS3",
            "{sample}_peaks.narrowPeak",
        ),
    params:
        out_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_peak_dir"],
            "MACS3",
            "{sample}",
        ),
        gs=config["genome_size"],
    singularity:
        "docker://aewebb/macs3:v3.0.1"
    resources:
        mem_mb=4000,
    threads: 4
    shell:
        "macs3 callpeak -t {input} -f BAMPE -g {params.gs} -n {params.out_prefix} --keep-dup all"
