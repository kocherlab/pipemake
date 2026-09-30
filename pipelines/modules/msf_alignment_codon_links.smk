ruleorder: link_tcs_score > link_trim_msa_clipkit > tcs_score > trim_msa_clipkit


use rule tcs_score as link_tcs_score with:
    input:
        "MSA/Codon/{sample}.fa",


use rule trim_msa_clipkit as link_trim_msa_clipkit with:
    input:
        "MSA/Codon/{sample}.fa",
