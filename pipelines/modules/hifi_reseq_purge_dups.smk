config = {
    "samples": ["bc2027"],
    "paths": {
        "workflow_prefix": "TEST",
        "reseq_fastq_dir": "fastq",
        "reseq_assembled_dir": "assembled",
    },
}


rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_assembled_dir"],
                "purge_dups",
                "{sample}.json",
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
        "purge_dups/purge_dups_v1.2.6.sif"
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
    run:
        import json

        with open(input[0], "r") as f:
            data = json.load(f)
        data["out_dir"] = params.out_dir
        with open(output[0], "w") as f:
            json.dump(data, f, indent=2)
