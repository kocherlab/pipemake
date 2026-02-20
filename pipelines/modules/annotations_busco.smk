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
        output_dir="Annotations/BUSCO",
        busco_db=config["busco_database"],
    singularity:
        "docker://ezlabgva/busco:v5.7.1_cv1"
    resources:
        mem_mb=12000,
    threads: 12
    shell:
        r"""
        busco -i {input} -o {params.output_dir} -l {params.busco_db} -m proteins -c {threads} --download_path {params.download_dir} -f &&
        summary_prefix=$(ls {params.output_dir} | grep -o 'short_summary.*BUSCO' | uniq) &&
        cp {params.output_dir}/${{summary_prefix}}.txt {params.output_dir}/summary.txt
        """
