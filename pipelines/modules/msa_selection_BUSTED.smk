rule all:
    input:
        expand(
            "Selection/BUSTED/{sample}.{label}.json",
            sample=config["samples"],
            label=config["busted_labels"],
        ),


rule label_tree_busted:
    input:
        gene_tree="iqtree/{sample}.fa.contree",
        species_tree="Selection/species_tree.tre",
    output:
        "Selection/BUSTED/{sample}.fa.{label}.tre",
    log:
        "logs/BUSTED/{sample}.{label}.log",
    params:
        output_dir=subpath(output[0], parent=True),
    singularity:
        "docker://aewebb/pipemake_utils:v1.4.3"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "label-gene-tree --species-tree {input.species_tree} --gene-tree {input.gene_tree} --analysis-label-str '{wildcards.label}:{wildcards.label}' --label-descendants --output-dir {params.output_dir} --output-naming-style suffix &> {log}"


rule busted_selection_hyphy:
    input:
        trimmed_msa="MSA/Trimmed/{sample}.fa",
        tree="Selection/BUSTED/{sample}.fa.{label}.tre",
    output:
        "Selection/BUSTED/{sample}.{label}.json",
    log:
        "logs/BUSTED/{sample}.{label}.log",
    singularity:
        "docker://aewebb/hyphy:v2.5.97_20260105"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "hyphy /hyphy/hyphy-analyses/BUSTED-MH/BUSTED-MH.bf --alignment {input.trimmed_msa} --tree {input.tree} --branches {wildcards.label} --output {output}"
