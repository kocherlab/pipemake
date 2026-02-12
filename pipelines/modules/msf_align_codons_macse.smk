rule all:
    input:
        expand("MSA/AA/{sample}.fa", sample=config["samples"]),
        expand("MSA/Codon/{sample}.fa", sample=config["samples"]),


rule msf_align_macse:
    input:
        "MSF/CDS/{sample}.fa",
    output:
        aa_msa="MSA/MACSE/{sample}_AA.macse_fmt.fa",
        nt_msa="MSA/MACSE/{sample}_NT.macse_fmt.fa",
    log:
        "logs/MACSE/{sample}.msf_align_macse.log",
    singularity:
        "docker://ranwez/omm_macse:v10.02"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "macse -prog alignSequences -seq {input} -out_NT {output.nt_msa} -out_AA {output.aa_msa}"


rule msa_export_macse:
    input:
        "MSA/MACSE/{sample}_NT.macse_fmt.fa",
    output:
        aa_msa="MSA/AA/{sample}.fa",
        nt_msa="MSA/Codon/{sample}.fa",
    log:
        "logs/MACSE/{sample}.msa_export_macse.log",
    params:
        codon_final_stop=config["macse_params"]["codon_final_stop"],
        codon_external_fs=config["macse_params"]["codon_external_fs"],
        codon_internal_fs=config["macse_params"]["codon_internal_fs"],
        char_remaining_fs=config["macse_params"]["char_remaining_fs"],
    singularity:
        "docker://ranwez/omm_macse:v10.02"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "macse -prog exportAlignment -align {input} -codonForFinalStop {params.codon_final_stop} -codonForExternalFS {params.codon_external_fs} -codonForInternalFS {params.codon_internal_fs}  -charForRemainingFS {params.char_remaining_fs} -out_NT {output.nt_msa} -out_AA {output.aa_msa}"
