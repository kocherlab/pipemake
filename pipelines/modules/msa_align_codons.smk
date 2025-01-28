rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["msa_codon_aligned_dir"],
                "{sample}.fasta",
            ),
            sample=config["samples"],
        ),


rule create_codon_aligned_msa:
    input:
        aligned_aa=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["msa_aligned_dir"],
            "MAFFT",
            "{sample}.fasta",
        ),
        unaligned_cds=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["msf_untranslated_dir"],
            "{sample}.fasta",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["msa_codon_aligned_dir"],
            "{sample}.fasta",
        ),
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.1"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "codon-alignment --aligned-aa {input.aligned_aa} --unaligned-cds {input.unaligned_cds} --output {output}"
