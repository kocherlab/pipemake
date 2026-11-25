rule all:
    input:
        "Annotations/BUSCO/summary.txt",


rule annotations_busco:
    input:
        f"Annotations/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
    output:
        "Annotations/BUSCO/summary.txt",
    log:
        "logs/BUSCO/annotations_busco.log",
    params:
        download_dir="Downloads/BUSCO",
        output_dir=subpath(output[0], parent=True),
        busco_db=config["busco_database"],
    singularity:
        "docker://ezlabgva/busco:v6.1.0_cv3"
    resources:
        mem_mb=12000,
    threads: 12
    shell:
        """
        busco -i {input} -o {params.output_dir} -l {params.busco_db} -m proteins -c {threads} --download_path {params.download_dir} -f &> {log}
        summary_prefix=$(ls {params.output_dir} | grep -o 'short_summary.*BUSCO' | uniq)
        cp {params.output_dir}/${{summary_prefix}}.txt {params.output_dir}/summary.txt
        """
