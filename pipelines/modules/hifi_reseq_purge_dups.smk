rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_assembled_dir"],
                "purge_dups",
                "{sample}",
                "seqs",
                "{sample}.p_ctg.purged.fa",
            ),
            sample=config["samples"],
        ),


rule create_fastq_list:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_assembled_dir"],
                "purge_dups",
                "{sample}.list",
            ),
        ),
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "echo {input} > {output}"


rule build_config:
    input:
        assembled_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "hifiasm",
            "{sample}.p_ctg.fa",
        ),
        fastq_list=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "purge_dups",
            "{sample}.list",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_assembled_dir"],
                "purge_dups",
                "{sample}.tmp.json",
            ),
        ),
    params:
        output_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "purge_dups",
            "{sample}_tmp",
        ),
    singularity:
        "docker://aewebb/purge_dups:v1.2.6"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "pd_config.py {input.assembled_fasta} {input.fastq_list} -n {output} -l {params.output_dir}"


rule update_json:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "purge_dups",
            "{sample}.tmp.json",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "purge_dups",
            "{sample}.json",
        ),
    params:
        out_dir=os.path.abspath(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_assembled_dir"],
                "purge_dups",
                "{sample}",
            ),
        ),
    resources:
        mem_mb=2000,
    threads: 1
    run:
        import json

        with open(input[0], "r") as f:
            data = json.load(f)
        data["out_dir"] = params.out_dir
        with open(output[0], "w") as f:
            json.dump(data, f, indent=2)


rule hifi_assembly_purge_dups:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "purge_dups",
            "{sample}.json",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "purge_dups",
            "{sample}",
            "seqs",
            "{sample}.p_ctg.purged.fa",
        ),
    params:
        species=config["species"],
    singularity:
        "docker://aewebb/purge_dups:v1.2.6"
    resources:
        mem_mb=30000,
    threads: 12
    shell:
        "run_purge_dups.py {input} /opt/conda/envs/purge_dups/bin {params.species} -p bash"


rule collect_purged_fasta:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "purge_dups",
            "{sample}",
            "seqs",
            "{sample}.p_ctg.purged.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            f"{{sample}}_{config['species']}_{config['assembly_version']}.fa",
        ),
    shell:
        "cp {input} {output}"
