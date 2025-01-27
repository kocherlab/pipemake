rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_aligned_bam_dir"],
                "{sample}.Aligned.bam",
            ),
            sample=config["samples"],
        ),


rule longread_reseq_align_minimap2:
    input:
        fasta_file=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
        reseq_fastqs=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_fastq_dir"],
            "{sample}_R1.fastq.gz",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_aligned_bam_dir"],
                "{sample}.Aligned.bam",
            ),
        ),
    singularity:
        "docker://aewebb/minimap2:v2.28"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "minimap2 -t {threads} -ax map-hifi -uf {input.fasta_file} {input.isoseq_fastqs} | samtools view --threads {threads} -bh -o {output}"
