rule all:
    input:
        "wASTRAL/species_tree.tre",


rule cat_gene_trees:
    input:
        expand("iqtree/{sample}.fa.contree", sample=config["samples"]),
    output:
        "wASTRAL/gene_trees.tre",
    localrule: True
    shell:
        "cat {input} > {output}"


rule infer_species_tree_wASTRAL:
    input:
        gene_trees="wASTRAL/gene_trees.tre",
        id_mappings="wASTRAL/id_mappings.txt",
    output:
        "wASTRAL/species_tree.tre",
    log:
        "logs/wASTRAL/species_tree.log",
    singularity:
        "docker://aewebb/aster:20260714"
    resources:
        mem_mb=32000,
    threads: 4
    shell:
        "wastral --input {input.gene_trees} --mapping {input.id_mappings} -t {threads} -output {output} 2> {log}"
