use rule trim_msa_clipkit as link_trim_msa_clipkit with:
    input:
        "MSA/NT/{sample}.fa",


use rule reconstruct_tree_iqtree as link_reconstruct_tree_iqtree with:
    input:
        "MSA/Trimmed/{sample}.fa",
