rule all:
    input:
        expand("MSA/NT/{sample}.fa", sample=config["samples"]),


use rule msf_align_muscle as link_msf_align_cds_muscle with:
    input:
        "MSF/NT/{sample}.fa",


rule msa_link_muscle_cds:
    input:
        "MSA/MUSCLE/{sample}.fa",
    output:
        "MSA/NT/{sample}.fa",
    localrule: True
    run:
        import os

        os.symlink(os.path.abspath(input[0]), output[0])
