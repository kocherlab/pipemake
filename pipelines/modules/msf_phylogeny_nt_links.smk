ruleorder: link_trim_msa_clipkit > link_create_iqtree_msa > trim_msa_clipkit > create_iqtree_msa


use rule trim_msa_clipkit as link_trim_msa_clipkit with:
    input:
        "MSA/NT/{sample}.fa",


use rule create_iqtree_msa as link_create_iqtree_msa with:
    input:
        "MSA/Trimmed/{sample}.fa",
