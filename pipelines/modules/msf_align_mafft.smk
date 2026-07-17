rule all:
    input:
        expand("MSA/MAFFT/{sample}.fa", sample=config["samples"]),


rule msf_align_mafft:
    input:
        "MSF/{sample}.fa",
    output:
        "MSA/MAFFT/{sample}.fa",
    params:
        run_mode='' if not config['mafft_params']['run_mode'] else f"-{config['mafft_params']['run_mode']}",
        op='' if not config['mafft_params']['op'] else f"--op {config['mafft_params']['op']}",
    log:
        "logs/MAFFT/{sample}.msf_align_mafft.log",
    singularity:
        "docker://aewebb/mafft:v7.525"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "mafft{params.run_mode} {params.op} {input} > {output} 2> {log}"
