rule all:
    input:
        expand(
            "ATAC_seq/BAM/Deduplicated/{sample}.deduplicated.bam",
            sample=config["samples"],
        ),


rule dedup_bam_sambamba:
    input:
        "ATAC_seq/BAM/Sorted/{sample}.sortedByCoord.bam",
    output:
        "ATAC_seq/BAM/Deduplicated/{sample}.deduplicated.bam",
    singularity:
        "docker://aewebb/sambamba:v1.0.1"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "sambamba markdup -t {threads} {input} {output}"
