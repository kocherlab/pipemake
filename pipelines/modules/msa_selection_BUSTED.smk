rule all:
    input:
        expand("Selection/BUSTED/{sample}.json", sample=config["samples"]),


rule busted_selection_hyphy:
    input:
        trimmed_msa="MSA/Trimmed/{sample}.fa",
        tree="iqtree/{sample}.fa.contree",
    output:
        "Selection/BUSTED/{sample}.json"
    singularity:
        "docker://aewebb/hyphy:v2.5.97_20260105"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "hyphy /hyphy/hyphy-analyses/BUSTED-MH/BUSTED-MH.bf --alignment {input.trimmed_msa} --tree {input.tree} --branches All --output {output}"
        
