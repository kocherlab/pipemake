rule all:
    input:
        expand("HiFi/Assembly/{sample}.purged.fa", sample=config["samples"]),


rule create_fastq_list:
    input:
        "HiFi/Assembly/{sample}.fa",
    output:
        temp("HiFi/Assembly/purge_dups/{sample}/file_input.list"),
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "echo {input} > {output}"


rule build_config:
    input:
        assembled_fasta="HiFi/Assembly/{sample}.fa",
        fastq_list="HiFi/Assembly/purge_dups/{sample}/file_input.list",
    output:
        json=temp("HiFi/Assembly/purge_dups/{sample}/.tmp.json"),
        tmp_dir=directory("HiFi/Assembly/purge_dups/{sample}/.tmp"),
    singularity:
        "docker://aewebb/purge_dups:v1.2.6"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        """
        pd_config.py {input.assembled_fasta} {input.fastq_list} -n {output.json} -l {output.tmp_dir}
        sleep 30
        """


rule update_json:
    input:
        "HiFi/Assembly/purge_dups/{sample}/.tmp.json",
    output:
        "HiFi/Assembly/purge_dups/{sample}/config.json",
    params:
        out_dir=subpath(output[0], parent=True),
        busco_db=config["busco_database"],
    resources:
        mem_mb=2000,
    threads: 1
    run:
        import json

        with open(input[0], "r") as f:
            data = json.load(f)
        data["out_dir"] = os.path.abspath(params.out_dir)
        data["busco"]["lineage"] = params.busco_db
        with open(output[0], "w") as f:
            json.dump(data, f, indent=2)


rule hifi_assembly_purge_dups:
    input:
        "HiFi/Assembly/purge_dups/{sample}/config.json",
    output:
        "HiFi/Assembly/purge_dups/{sample}/seqs/{sample}.purged.fa",
    params:
        sample="{sample}",
    singularity:
        "docker://aewebb/purge_dups:v1.2.6"
    resources:
        mem_mb=30000,
    threads: 12
    shell:
        "run_purge_dups.py {input} /opt/conda/envs/purge_dups/bin {params.sample} -p bash"


rule collect_purged_fasta:
    input:
        "HiFi/Assembly/purge_dups/{sample}/seqs/{sample}.purged.fa",
    output:
        "HiFi/Assembly/{sample}.purged.fa",
    params:
        output_dir=subpath(input[0], parent=True) + "/.tmp",
    localrule: True
    shell:
        """
        cp {input} {output}
        rm -rf {params.output_dir}
        """
