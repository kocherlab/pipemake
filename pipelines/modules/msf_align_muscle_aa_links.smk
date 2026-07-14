rule all:
    input:
        expand("MSA/AA/{sample}.fa", sample=config["samples"]),


use rule msf_align_muscle as link_msf_align_aa_muscle with:
    input:
        "MSF/AA/{sample}.fa",


rule msa_link_muscle_aa:
    input:
        "MSA/MUSCLE/{sample}.fa",
    output:
        "MSA/AA/{sample}.fa",
    localrule: True
    run:
        import os

        os.symlink(os.path.abspath(input[0]), output[0])
