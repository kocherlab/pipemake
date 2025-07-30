rule all:
    input:
        f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa",


rule hifi_wo_hic_reseq_assemble_hifiasm:
    input:
        f"HiFi/FASTQ/{config['species']}.fq.gz",
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "hifiasm",
            "{sample}.a_ctg.gfa",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["reseq_assembled_dir"],
            "hifiasm",
            "{sample}.p_ctg.gfa",
        ),
#        f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.p_ctg.gfa",
    params:
        output_prefix=f"Assembly/hifiasm/{config['species']}",
    singularity:
        "docker://aewebb/hifiasm:v0.24.0-r702"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        "hifiasm --primary -t {threads} -o {params.output_prefix} {input}"


rule hifi_cvt_hifiasm_gfa_to_fasta:
    input:
        f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.p_ctg.gfa",
    output:
        f"Assembly/hifiasm/{config['species']}_{config['assembly_version']}.fa",
    shell:
        """
        awk '/^S/{{print ">"$2"\\n"$3}}' {input} | fold > {output}
        sleep 30
        """
