rule all:
    input:
        expand(f"Selection/BUSTED/{{sample}}.{config['busted_label']}.json", sample=config["samples"]),


rule label_tree_busted:
    input:
        gene_tree="iqtree/{sample}.fa.contree",
        species_tree="Selection/species_tree.tre",
    output:
        f"Selection/BUSTED/{{sample}}.fa.{config['busted_label']}.tre"
    log:
        f"logs/BUSTED/{{sample}}.label.log",
    params:
        output_dir=subpath(output[0], parent=True),
        label=config['busted_label'],
        node=config["busted_node"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.4.2"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "label-gene-tree --species-tree {input.species_tree} --gene-tree {input.gene_tree} --analysis-label-str '{params.label}:{params.node}' --label-descendants --output-dir {params.output_dir} --output-naming-style suffix &> {log}"

rule busted_selection_hyphy:
    input:
        trimmed_msa="MSA/Trimmed/{sample}.fa",
        tree=f"Selection/BUSTED/{{sample}}.fa.{config['busted_label']}.tre",
    output:
        f"Selection/BUSTED/{{sample}}.{config['busted_label']}.json",
    log:
        f"logs/BUSTED/{{sample}}.log",
    params:
        label=config['busted_label'],
    singularity:
        "docker://aewebb/hyphy:v2.5.97_20260105"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "hyphy /hyphy/hyphy-analyses/BUSTED-MH/BUSTED-MH.bf --alignment {input.trimmed_msa} --tree {input.tree} --branches {params.label} --output {output}"