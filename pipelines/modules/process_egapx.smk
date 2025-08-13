rule all:
    input:
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gtf",
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa",
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
        f"Processed/{config['species']}_genome_{config['assembly_version']}.fasta.gz",
        f"Processed/{config['species']}_genome_{config['assembly_version']}.fasta.gz.fai",
        f"Processed/{config['species']}_genome_{config['assembly_version']}.fasta.gz.gzi",
        f"Processed/{config['species']}_genome_{config['assembly_version']}.assembly.stats",


rule process_egapx_gtf:
    input:
        config["gtf_input"],
    output:
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gtf",
        temp(
            f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.no_utrs.gff"
        ),
    params:
        out_prefix=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}",
        species=config["species"],
        assembly_version=config["assembly_version"],
        annotation_version=config["annotation_version"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.6"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        """
        process-ncbi-annotations --gtf {input} --species-tag {params.species}_{params.assembly_version}.{params.annotation_version} --out-prefix {params.out_prefix}
        mv {params.out_prefix}.gff {params.out_prefix}.no_utrs.gff
        """


rule add_utrs_to_gff:
    input:
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.no_utrs.gff",
    output:
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
    singularity:
        "docker://aewebb/ncbi-genome-tools:20250625"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "add_utrs_to_gff.py {input} > {output}"


rule assembly_stats_bbmap:
    input:
        config["assembly_fasta"],
    output:
        f"Processed/{config['species']}_genome_{config['assembly_version']}.assembly.stats",
    singularity:
        "docker://biocontainers/bbmap:39.26--he5f24ec_0"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "stats.sh in={input} > {output}"


rule gff_to_transcripts:
    input:
        gff_file=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        assembly_fasta=config["assembly_fasta"],
    output:
        temp(
            f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa.tmp"
        ),
    singularity:
        "docker://aewebb/gffread:v0.12.7"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gffread -x {output} -g {input.assembly_fasta} {input.gff_file}"


rule update_transcripts:
    input:
        gff_file=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        transcript_fasta=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa.tmp",
    output:
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa",
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.6"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "update-fasta --gff-file {input.gff_file} --fasta-file {input.transcript_fasta} --out-file {output} --seq-type nucleotide"


rule gff_to_proteins:
    input:
        gff_file=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        assembly_fasta=config["assembly_fasta"],
    output:
        temp(
            f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa.tmp"
        ),
    singularity:
        "docker://aewebb/gffread:v0.12.7"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gffread -y {output} -g {input.assembly_fasta} {input.gff_file}"


rule update_proteins:
    input:
        gff_file=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        protein_fasta=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa.tmp",
    output:
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa",
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.6"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "update-fasta --gff-file {input.gff_file} --fasta-file {input.protein_fasta} --out-file {output} --seq-type protein"


rule gzip_assembly:
    input:
        config["assembly_fasta"],
    output:
        f"Processed/{config['species']}_genome_{config['assembly_version']}.fasta.gz",
    singularity:
        "docker://aewebb/tabix:v1.11"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "bgzip -c {input} > {output}"


rule index_assembly:
    input:
        f"Processed/{config['species']}_genome_{config['assembly_version']}.fasta.gz",
    output:
        f"Processed/{config['species']}_genome_{config['assembly_version']}.fasta.gz.fai",
        f"Processed/{config['species']}_genome_{config['assembly_version']}.fasta.gz.gzi",
    singularity:
        "docker://aewebb/samtools:v1.20"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "samtools faidx {input}"
