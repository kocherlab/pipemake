rule all:
    input:
        expand("MSA/NT/{sample}.fa", sample=config["samples"]),


use rule msf_align_mafft as link_msf_align_nt_mafft with:
    input:
        "MSF/NT/{sample}.fa",


rule msa_link_mafft_nt:
    input:
        "MSA/MAFFT/{sample}.fa",
    output:
        "MSA/NT/{sample}.fa",
    localrule: True
    run:
        import os

        os.symlink(os.path.abspath(input[0]), output[0])
