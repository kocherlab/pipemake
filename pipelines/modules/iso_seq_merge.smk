rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["isoseq_merged_bam_dir"],
            f"{config['species']}_{config['assembly_version']}.bam",
        ),


rule merge_bam:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["isoseq_sorted_bam_dir"],
                "{sample}.sortedByCoord.bam",
            ),
            sample=config["samples"],
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["isoseq_merged_bam_dir"], 
            f"{config['species']}_{config['assembly_version']}.bam",
        ),
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "samtools merge -@ {threads} -r {output} {input}"

        singularity exec -B ${PWD}:${PWD} braker3_lr.sif braker.pl --genome=genome.fa --prot_seq=protein_db.fa –-bam=isoseq.bam --threads=${T}
