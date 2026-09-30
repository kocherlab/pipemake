ruleorder: link_tcs_score > link_trim_msa_clipkit > link_create_iqtree_msa > tcs_score > trim_msa_clipkit > create_iqtree_msa


use rule tcs_score as link_tcs_score with:
    input:
        "MSA/AA/{sample}.fa",


use rule trim_msa_clipkit as link_trim_msa_clipkit with:
    input:
        "MSA/AA/{sample}.fa",


use rule create_iqtree_msa as link_create_iqtree_msa with:
    input:
        "MSA/Trimmed/{sample}.fa",
