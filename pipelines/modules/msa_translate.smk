rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["msf_unaligned_dir"],
                "{sample}.fasta",
            ),
            sample=config["samples"],
        ),


rule translate_msf:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["msf_untranslated_dir"],
            "{sample}.fasta",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["msf_unaligned_dir"],
            "{sample}.fasta",
        ),
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.0"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "translate-seq-file --input {input} --output {output}"
