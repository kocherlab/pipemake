rule all:
    input:
        expand("MSA/Codon/{sample}.fasta", sample=config["samples"]),


rule create_codon_aligned_msa:
    input:
        aligned_aa="MSA/AA/{sample}.fasta",
        unaligned_cds="MSF/CDS/{sample}.fasta",
    output:
        "MSA/Codon/{sample}.fasta",
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "codon-alignment --aligned-aa {input.aligned_aa} --unaligned-cds {input.unaligned_cds} --output {output}"
