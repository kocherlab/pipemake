rule all:
    input:
        expand("Selection/absrel/{sample}.json", sample=config["samples"]),


rule absrel_selection_hyphy:
    input:
        trimmed_msa="MSA/Trimmed/{sample}.fa",
        tree="iqtree/{sample}.fa.contree",
    output:
        "Selection/absrel/{sample}.json"
    singularity:
        "docker://aewebb/hyphy:v2.5.97_20260105"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "hyphy absrel --alignment {input.trimmed_msa} --tree {input.tree} --branches All --multiple-hits Double+Triple --output {output}"
