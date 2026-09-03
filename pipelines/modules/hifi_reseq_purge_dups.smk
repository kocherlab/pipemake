rule all:
    input:
        f"Assembly/purge_dups/{config['species']}_{config['assembly_version']}.fa",


rule create_fastq_list:
    input:
        f"HiFi/FASTQ/{config['species']}.fq.gz",
    output:
        temp(f"Assembly/purge_dups/{config['species']}.list"),
    localrule: True
    shell:
        "echo {input} > {output}"


rule build_config:
    input:
        assembled_fasta=f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa",
        fastq_list=f"Assembly/purge_dups/{config['species']}.list",
    output:
        temp(f"Assembly/purge_dups/{config['species']}.tmp.json"),
    log:
        f"logs/purge_dups/{config['species']}.log",
    params:
        output_dir=subpath(output[0], strip_suffix=".tmp.json") + "_tmp",
    singularity:
        "docker://aewebb/purge_dups:v1.2.6"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        """
        pd_config.py {input.assembled_fasta} {input.fastq_list} -n {output} -l {params.output_dir} &> {log}
        sleep 30
        """


rule update_json:
    input:
        f"Assembly/purge_dups/{config['species']}.tmp.json",
    output:
        f"Assembly/purge_dups/{config['species']}.json",
    params:
        out_dir=os.path.abspath(subpath(output[0], strip_suffix=".json")),
        busco_db=config["busco_database"],
    localrule: True
    run:
        import json

        with open(input[0], "r") as f:
            data = json.load(f)
        data["out_dir"] = params.out_dir
        data["busco"]["lineage"] = params.busco_db
        with open(output[0], "w") as f:
            json.dump(data, f, indent=2)


rule hifi_assembly_purge_dups:
    input:
        f"Assembly/purge_dups/{config['species']}.json",
    output:
        f"Assembly/purge_dups/{config['species']}/seqs/{config['species']}_{config['assembly_version']}.purged.fa",
    log:
        f"logs/purge_dups/{config['species']}.log",
    params:
        species=config["species"],
    singularity:
        "docker://aewebb/purge_dups:v1.2.6"
    resources:
        mem_mb=30000,
    threads: 12
    shell:
        "run_purge_dups.py {input} /opt/conda/envs/purge_dups/bin {params.species} -p bash &> {log}"


rule collect_purged_fasta:
    input:
        f"Assembly/purge_dups/{config['species']}/seqs/{config['species']}_{config['assembly_version']}.purged.fa",
    output:
        f"Assembly/purge_dups/{config['species']}_{config['assembly_version']}.fa",
    params:
        output_dir=subpath(output[0], strip_suffix=f"_{config['assembly_version']}.fa") + "_tmp",
    localrule: True
    shell:
        """
        cp {input} {output}
        rm -rf {params.output_dir}
        """
