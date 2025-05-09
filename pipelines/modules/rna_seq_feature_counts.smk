rule all:
    input:
        f"RNAseq/Counts/{config['species']}_{config['assembly_version']}.featurecounts.tsv",


ruleorder: feature_counts_pair_end > feature_counts_single_end


rule feature_counts_pair_end:
    input:
        gtf_file=f"Assembly/{config['species']}_{config['assembly_version']}.gtf",
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        r2_reads="RNAseq/FASTQ/{sample}_R2.fq.gz",
        sorted_bam="RNAseq/BAM/Sorted/{sample}.sortedByCoord.bam",
    output:
        "RNAseq/Counts/featureCounts/{sample}.featurecounts.txt",
    params:
        count_strand=config["count_strand"],
    singularity:
        "docker://aewebb/subread:v2.0.1"
    resources:
        mem_mb=4000,
    threads: 4
    shell:
        "featureCounts -s {params.count_strand} -p -T {threads} -a {input.gtf_file} -o {output} {input.sorted_bam}"


rule feature_counts_single_end:
    input:
        gtf_file=f"Assembly{config['species']}_{config['assembly_version']}.gtf",
        r1_reads="RNAseq/FASTQ/{sample}_R1.fq.gz",
        sorted_bam="RNAseq/BAM/Sorted/{sample}.sortedByCoord.bam",
    output:
        "RNAseq/Counts/featureCounts/{sample}.featurecounts.txt",
    params:
        count_strand=config["count_strand"],
    singularity:
        "docker://aewebb/subread:v2.0.1"
    resources:
        mem_mb=4000,
    threads: 4
    shell:
        "featureCounts -s {params.count_strand} -T {threads} -a {input.gtf_file} -o {output} {input.sorted_bam}"


rule feature_counts_report:
    input:
        expand(
            "RNAseq/Counts/featureCounts/{sample}.featurecounts.txt",
            sample=config["samples"],
        ),
    output:
        f"RNAseq/Counts/{config['species']}_{config['assembly_version']}.featurecounts.tsv",
    params:
        count_dir="RNAseq/Counts/featureCounts",
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "featureCounts-report {params.count_dir} {output}"
