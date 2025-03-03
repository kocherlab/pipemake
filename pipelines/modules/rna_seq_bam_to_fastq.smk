ruleorder: bam_to_pair_end_fastq > bam_to_single_end_fastq


rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["rnaseq_filtered_fastq_dir"],
                ".{sample}.chk",
            ),
            sample=config["samples"],
        ),


rule bam_to_single_end_fastq:
    input:
        bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_sorted_bam_dir"],
            "{sample}.sortedByCoord.bam",
        ),
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    output:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_filtered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        cvt_chk=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["rnaseq_filtered_fastq_dir"],
                ".{sample}.chk",
            )
        ),
    singularity:
        "docker://aewebb/bedtools:v2.31.1"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "bedtools bamtofastq -i {input.bam} -fq {output.r1_reads} && touch {output.cvt_chk}"


rule bam_to_pair_end_fastq:
    input:
        bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_sorted_bam_dir"],
            "{sample}.sortedByCoord.bam",
        ),
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        r2_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_fastq_dir"],
            "{sample}_R2.fq.gz",
        ),
    output:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_filtered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
        r2_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq_filtered_fastq_dir"],
            "{sample}_R2.fq.gz",
        ),
        cvt_chk=temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["rnaseq_filtered_fastq_dir"],
                ".{sample}.chk",
            )
        ),
    singularity:
        "docker://aewebb/bedtools:v2.31.1"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "bedtools bamtofastq -i {input.bam} -fq {output.r1_reads} -fq2 {output.r2_reads} && touch {output.cvt_chk}"
