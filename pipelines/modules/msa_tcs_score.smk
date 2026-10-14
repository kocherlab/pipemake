config-{
  "samples": ["test1"],
  "tcoffee_params": {
    "evaluate_mode": "fmcoffee"
  }
}

rule all:
    input:
        expand("Stats/MSA/{sample}.tcs_score.txt", sample=config["samples"]),
        expand("Stats/MSA/{sample}.score_ascii", sample=config["samples"]),
        expand("Stats/MSA/{sample}.html", sample=config["samples"]),


rule tcs_score:
    input:
        "MSA/{sample}.fa",
    output:
        tcs_score="Stats/MSA/{sample}.tcs_score.txt",
        ascii_results="Stats/MSA/{sample}.score_ascii",
        html_results="Stats/MSA/{sample}.html",
    log:
        "logs/t_coffee/{sample}.tcs.log",
    params:
        out_prefix=subpath(output.ascii_results, strip_suffix=".score_ascii"),
        evaluate_mode='-mode=fmcoffee' if 'fmcoffee' in config['tcoffee_params']['evaluate_mode'] else '',
    singularity:
        "docker://aewebb/tcoffee:v13.41.0.28bdc39"
    resources:
        mem_mb=16000,
    threads: 1
    shell:  
        """
        t_coffee -infile={input} -evaluate {params.evaluate_mode} -outfile={params.out_prefix} -no_warning -output score_ascii &> {log}
        sed -n 's/^SCORE=//p' {output.ascii_results} | awk '{{printf \"%.2f\\n\", $1/100}}' > {output.tcs_score}
        rm {params.out_prefix}
        """
