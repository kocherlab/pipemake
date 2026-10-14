rule all:
    input:
        expand("Stats/MSA/{sample}.tcs_score.txt", sample=config["samples"]),
        expand("Stats/MSA/{sample}.score_ascii", sample=config["samples"]),


rule tcs_score:
    input:
        "MSA/{sample}.fa",
    output:
        tcs_score="Stats/MSA/{sample}.tcs_score.txt",
        ascii_results="Stats/MSA/{sample}.score_ascii",
    log:
        "logs/t_coffee/{sample}.tcs.log",
    singularity:
        "docker://aewebb/tcoffee:v13.41.0.28bdc39"
    resources:
        mem_mb=16000,
    threads: 1
    shell:  
        """
        t_coffee -infile {input} -evaluate -mode fmcoffee -outfile {output.ascii_results} -no_warning -output score_ascii &> {log}
        sed -n 's/^SCORE=//p' {output.ascii_results} | awk '{{printf \"%.2f\\n\", $1/100}}' > {output.tcs_score}
        rm {params.out_prefix}
        """
