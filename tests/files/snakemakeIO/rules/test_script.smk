rule script_call:
    input:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    output:
        r1_reads=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["filtered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    resources:
        mem_mb=config["resources"]["script_call"]["mem_mb"],
    threads: config["resources"]["script_call"]["threads"]
    script:
        "../scripts/test.py"
