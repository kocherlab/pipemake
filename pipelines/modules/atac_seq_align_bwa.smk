ruleorder: bwa_mem_pair_end_atac_seq > bwa_mem_single_end_atac_seq


rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["atacseq_aligned_bam_dir"],
                "{sample}.Aligned.bam",
            ),
            sample=config["samples"],
        ),


rule bwa_index_atac_seq:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "BWA",
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "BWA",
            f"{config['species']}_{config['assembly_version']}.fa.amb",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "BWA",
            f"{config['species']}_{config['assembly_version']}.fa.ann",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "BWA",
            f"{config['species']}_{config['assembly_version']}.fa.bwt",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "BWA",
            f"{config['species']}_{config['assembly_version']}.fa.pac",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "BWA",
            f"{config['species']}_{config['assembly_version']}.fa.sa",
        ),
    params:
        index_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "BWA",
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
    singularity:
        "docker://aewebb/bwa:v0.7.18"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "cp {input} {params.index_fasta} && "
        "bwa index {params.index_fasta}"


rule bwa_mem_single_end_atac_seq:
    input:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_fastq_dir"],
            "{sample}-R1.fastq.gz",
        ),
        assembly=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
        index_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "BWA",
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_aligned_bam_dir"],
            "{sample}.Aligned.bam",
        ),
    singularity:
        "docker://aewebb/bwa:v0.7.18"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "bwa mem -t {threads} {input.index_fasta} {input.r1_reads} | samtools view --threads {threads} -bh -o {output}"


rule bwa_mem_pair_end_atac_seq:
    input:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_fastq_dir"],
            "{sample}-R1.fastq.gz",
        ),
        r2_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_fastq_dir"],
            "{sample}-R2.fastq.gz",
        ),
        index_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["index_dir"],
            "BWA",
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["atacseq_aligned_bam_dir"],
            "{sample}.Aligned.bam",
        ),
    singularity:
        "docker://aewebb/bwa:v0.7.18"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "bwa mem -t {threads} {input.index_fasta} {input.r1_reads} {input.r2_reads} | samtools view --threads {threads} -bh -o {output}"