rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.gff3",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.codingseq",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.aa",
        ),


rule annotate_braker3:
    input:
        masked_assembly=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa.masked",
        ),
        merged_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["rnaseq-merged-bam-dir"],
            f"{config['species']}_{config['assembly_version']}.bam",
        ),
        protein_hints=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["homology_dir"],
            "ProteinHints.fa",
        ),
        augustus_check=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["downloads_dir"],
            "augustus",
            f".config.chk",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.gff3",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.codingseq",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
            "braker.aa",
        ),
    params:
        annotations_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["annotations_dir"],
            "BRAKER3",
        ),
        augustus_config=os.path.abspath(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["downloads_dir"],
                "augustus",
                "config",
            )
        ),
    singularity:
        "docker://teambraker/braker3:v3.0.7.6"
    resources:
        mem_mb=32000,
    threads: 20
    shell:
        "braker.pl --genome {input.masked_assembly} --prot_seq {input.protein_hints} --bam {input.merged_bam} -gff3 --softmasking --threads {threads} --workingdir {params.annotations_dir} --AUGUSTUS_CONFIG_PATH {params.augustus_config}"
