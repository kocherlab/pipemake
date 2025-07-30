rule all:
    input:
        "Annotations/compleasm/summary.txt",


rule annotations_compleasm:
    input:
        transcript_fasta=f"Annotations/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa",
        compleasm_library=f"Downloads/compleasm/{config['busco_database']}.done",
    output:
        "Annotations/compleasm/summary.txt",
    params:
        output_dir="Annotations/compleasm",
        download_dir="Downloads/compleasm",
        busco_db=config["busco_database"],
    singularity:
        "docker://huangnengcsu/compleasm:v0.2.6"
    resources:
        mem_mb=12000,
    threads: 12
    shell:
        "compleasm run -t{threads} -L {params.download_dir} -l {params.busco_db} -a {input.transcript_fasta} -o {params.output_dir}"
