rule all:
    input:
        expand("Selection/RELAX/{sample}.json", sample=config["samples"]),


rule relax_selection_hyphy:
    input:
        trimmed_msa="MSA/Trimmed/{sample}.fa",
        tree="iqtree/{sample}.fa.contree",
    output:
        "Selection/RELAX/{sample}.json"
    singularity:
        "docker://aewebb/hyphy:v2.5.97_20260105"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "hyphy relax --alignment {input.trimmed_msa} --tree {input.tree} --multiple-hits Double+Triple --output {output}"
        
