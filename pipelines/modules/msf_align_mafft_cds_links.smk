rule all:
    input:
        expand("MSA/CDS/{sample}.fa", sample=config["samples"]),


use rule msf_align_mafft as link_msf_align_cds_mafft with:
    input:
        "MSF/CDS/{sample}.fa",


rule msa_link_mafft_cds:
    input:
        "MSA/MAFFT/{sample}.fa",
    output:
        "MSA/CDS/{sample}.fa",
    localrule: True
    shell:
        "ln -s {input} {output}"
