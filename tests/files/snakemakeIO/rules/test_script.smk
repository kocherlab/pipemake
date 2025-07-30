rule script_call:
    input:
        r1_reads=os.path.join(
            config["paths"]["unfiltered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    output:
        r1_reads=os.path.join(
            config["paths"]["filtered_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    resources:
        mem_mb=16000,
    threads: 4
    script:
        "../scripts/test.py"
