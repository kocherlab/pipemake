rule all:
    input:
        expand(
            "Assembly/hifiasm/{sample}.bp.{hap}.p_ctg.fa",
            sample=config["samples"],
            hap=["hap1", "hap2"],
        ),


rule hifi_wo_hic_reseq_assemble_phased_hifiasm:
    input:
        "HiFi/FASTQ/{sample}.fq.gz",
    output:
        "Assembly/hifiasm/{sample}.bp.hap1.p_ctg.gfa",
        "Assembly/hifiasm/{sample}.bp.hap2.p_ctg.gfa",
    params:
        output_prefix=subpath(output[0], strip_suffix=".bp.hap1.p_ctg.gfa"),
    singularity:
        "docker://aewebb/hifiasm:v0.24.0-r702"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        "hifiasm -t {threads} -o {params.output_prefix} {input}"


rule hifi_cvt_phased_hifiasm_gfa_to_fasta:
    input:
        "Assembly/hifiasm/{sample}.bp.{hap}.p_ctg.gfa",
    output:
        "Assembly/hifiasm/{sample}.bp.{hap}.p_ctg.fa",
    shell:
        """
        awk '/^S/{{print ">"$2"\\n"$3}}' {input} | fold > {output}
        sleep 30
        """
