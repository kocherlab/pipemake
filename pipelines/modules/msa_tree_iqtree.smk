rule all:
    input:
        expand("iqtree/{sample}.fa.contree", sample=config["samples"]),


rule create_iqtree_msa:
    input:
        "MSA/{sample}.fa",
    output:
        "iqtree/{sample}.fa",
    localrule: True
    run:
        import shutil

        shutil.copy(input[0], output[0])


ruleorder: reconstruct_tree_iqtree_with_outgroup > reconstruct_tree_iqtree


rule reconstruct_tree_iqtree:
    input:
        "iqtree/{sample}.fa",
    output:
        "iqtree/{sample}.fa.contree",
    log:
        "logs/iqtree/{sample}.log",
    params:
        msa=subpath(output[0], strip_suffix=".contree"),
    singularity:
        "docker://aewebb/iqtree:v3.0.1"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "iqtree -s {params.msa} -alrt 1000 -B 1000 -T {threads} &> {log}"


rule create_outgroup_file:
    input:
        tree_msa="iqtree/{sample}.fa",
        outgroup_file="Selection/outgroup_species.txt",
    output:
        "iqtree/{sample}.outgroup",
    singularity:
        "docker://aewebb/pipemake_utils:v1.4.3"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "confirm-ids-fasta --fasta {input.tree_msa} --ids-file {input.outgroup_file} --prefix-ids --first-only --output {output}"


rule reconstruct_tree_iqtree_with_outgroup:
    input:
        tree_msa="iqtree/{sample}.fa",
        outgroup="iqtree/{sample}.outgroup",
    output:
        "iqtree/{sample}.fa.contree",
    log:
        "logs/iqtree/{sample}.log",
    params:
        msa=subpath(output[0], strip_suffix=".contree"),
    singularity:
        "docker://aewebb/iqtree:v3.0.1"
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "iqtree -s {params.msa} -o {input.outgroup} -alrt 1000 -B 1000 -T {threads} &> {log}"
