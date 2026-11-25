rule all:
    input:
        "Assembly/compleasm/summary.txt",


rule assembly_compleasm:
    input:
        assembly_fasta=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
        compleasm_library=f"Downloads/compleasm/{config['busco_database']}.done",
    output:
        "Assembly/compleasm/summary.txt",
    log:
        "logs/compleasm/assembly_compleasm.log",
    params:
        output_dir=subpath(output[0], parent=True),
        download_dir=subpath(input.compleasm_library, parent=True),
        busco_db=config["busco_database"],
    singularity:
        "docker://aewebb/compleasm:v0.2.10"
    resources:
        mem_mb=12000,
    threads: 12
    shell:
        "compleasm run -t{threads} -L {params.download_dir} -l {params.busco_db} -a {input.assembly_fasta} -o {params.output_dir} &> {log}"
