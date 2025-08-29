rule all:
    input:
        expand("MSA/PRANK/{sample}.fa", sample=config["samples"]),


rule msf_align_codons_prank:
    input:
        "MSF/CDS/{sample}.fa",
    output:
        "MSA/PRANK/{sample}.fa",
    params:
        out_prefix="MSA/PRANK/{sample}",
    singularity:
        "docker://aewebb/prank:170703"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        """
        prank -codon -F -d={input} -o={params.out_prefix}
        mv {params.out_prefix}.best.fas {output}
        sleep 10
        """
