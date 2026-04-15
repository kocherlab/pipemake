rule all:
    input:
        expand(
            "reSEQ/PopGen/Admixture/log_{k}.out", k=range(1, int(config["max_k"]) + 1)
        ),


rule reseq_link_bed_files_for_admixture:
    input:
        bed_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bed",
        fam_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.fam",
    output:
        bed_file=temp(
            f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.admixture.bed"
        ),
        fam_file=temp(
            f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.admixture.fam"
        ),
    localrule: True
    shell:
        "ln -s {input.bed_file} {output.bed_file} && ln -s {input.fam_file} {output.fam_file}"


rule reseq_create_bim_file_for_admixture:
    input:
        f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.bim",
    output:
        temp(
            f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.admixture.bim"
        ),
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        awk 'BEGIN{{FS=\"[ \\t]+\"}}{{delim=match($0,/\\t/)?\"\\t\":\" \";if(!($1 in map)){{map[$1]=++i}};$1=map[$1];out=$1;for(j=2;j<=NF;j++)out=out delim $j;print out}}' {input} > {output}
        """


rule reseq_run_admixture:
    input:
        bed_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.admixture.bed",
        fam_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.admixture.fam",
        bim_file=f"reSEQ/PLINK/Filtered/{config['species']}_{config['assembly_version']}.filtered.admixture.bim",
    output:
        "reSEQ/PopGen/Admixture/log_{k}.out",
    resources:
        mem_mb=16000,
    threads: 4
    shell:
        "admixture --cv -j{threads} {input.bed_file} {wildcards.k} > reSEQ/PopGen/Admixture/log_{wildcards.k}.out"
