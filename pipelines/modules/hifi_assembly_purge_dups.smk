rule all:
    input:
        expand("HiFi/Assembly/{sample}.purged.fa", sample=config["samples"]),


rule create_fastq_list:
    input:
        "HiFi/Assembly/{sample}.fa",
    output:
        temp("HiFi/Assembly/purge_dups/{sample}.list"),
    localrule: True
    shell:
        "echo {input} > {output}"


rule build_config:
    input:
        assembled_fasta="HiFi/Assembly/{sample}.fa",
        fastq_list="HiFi/Assembly/purge_dups/{sample}.list",
    output:
        json=temp("HiFi/Assembly/purge_dups/{sample}.tmp.json"),
        fofn=temp("HiFi/Assembly/purge_dups/{sample}_tmp/pb.fofn"),
    params:
        out_dir=subpath(output.fofn, parent=True),
    singularity:
        "docker://aewebb/purge_dups:v1.2.6"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        """
        pd_config.py {input.assembled_fasta} {input.fastq_list} -n {output.json} -l {params.out_dir}
        sleep 30
        """


rule update_json:
    input:
        "HiFi/Assembly/purge_dups/{sample}.tmp.json",
    output:
        "HiFi/Assembly/purge_dups/{sample}.json",
    params:
        out_dir=subpath(output[0], parent=True),
        busco_db=config["busco_database"],
    resources:
        mem_mb=2000,
    threads: 1
    run:
        import os
        import json

        with open(input[0], "r") as f:
            data = json.load(f)
        data["out_dir"] = os.path.abspath(params.out_dir)
        data["busco"]["lineage"] = params.busco_db
        with open(output[0], "w") as f:
            json.dump(data, f, indent=2)


rule hifi_assembly_purge_dups:
    input:
        json="HiFi/Assembly/purge_dups/{sample}.json",
        fofn="HiFi/Assembly/purge_dups/{sample}_tmp/pb.fofn",
    output:
        "HiFi/Assembly/purge_dups/seqs/{sample}.purged.fa",
    params:
        tmp_dir=subpath(input.fofn, parent=True),
    singularity:
        "docker://aewebb/purge_dups:v1.2.6"
    resources:
        mem_mb=30000,
    threads: 12
    shell:
        """
        run_purge_dups.py {input.json} /opt/conda/envs/purge_dups/bin {wildcard.sample} -p bash
        rm -rf {params.tmp_dir}
        """


rule collect_purged_fasta:
    input:
        "HiFi/Assembly/purge_dups/seqs/{sample}.purged.fa",
    output:
        "HiFi/Assembly/{sample}.purged.fa",
    shell:
        "cp {input} {output}"
