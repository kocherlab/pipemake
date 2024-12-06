rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["msa_aligned_dir"],
                "MAFFT",
                "{sample}.fasta",
            ),
            sample=config["samples"],
        ),


rule align_msa_mafft:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["msf_unaligned_dir"],
            "{sample}.fasta",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["msa_aligned_dir"],
            "MAFFT",
            "{sample}.fasta",
        ),
    singularity:
        "docker://aewebb/mafft:v7.525"
    resources:
        mem_mb=8000,
    threads: 4
    shell:
        "mafft-linsi --thread {threads} {input} > {output}"
