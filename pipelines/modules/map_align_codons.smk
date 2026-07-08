rule all:
    input:
        expand("MSA/Codon/{sample}.fa", sample=config["samples"]),


rule create_codon_aligned_msa:
    input:
        aligned_aa="MSA/AA/{sample}.fa",
        unaligned_cds="MSF/CDS/{sample}.fa",
    output:
        "MSA/Codon/{sample}.fa",
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.8"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "codon-alignment --aligned-aa {input.aligned_aa} --unaligned-cds {input.unaligned_cds} --output {output}"
