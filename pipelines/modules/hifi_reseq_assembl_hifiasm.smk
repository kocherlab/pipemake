rule all:
    input:
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["reseq_assembled_dir"],
                "hifiasm",
                "{sample}.p_ctg.fa",
            ),
            sample=config["samples"],
        ),


rule hifi_wo_hic_reseq_assemble_hifiasm:
    input:
        reseq_fastq=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_fastq_dir"],
            "{sample}_R1.fq.gz",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "hifiasm",
            "{sample}.p_ctg.gfa",
        ),
    params:
        output_prefix=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "{sample}",
        ),
    singularity:
        "docker://aewebb/hifiasm:v0.24.0-r702"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        "hifiasm --primary -t {threads} -o {params.output_prefix} {input.reseq_fastq}"


rule hifi_cvt_hifiasm_gfa_to_fasta:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "hifiasm",
            "{sample}.p_ctg.gfa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "hifiasm",
            "{sample}.p_ctg.fa",
        ),
    shell:
        """
        awk '/^S/{{print ">"$2"\\n"$3}}' {input} | fold > {output}
        """
