rule all:
    input:
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gtf",
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa",
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa",

use rule assembly_compleasm as link_assembly_compleasm with:
    input:
        assembly_fasta=config["assembly_fasta"],

rule add_kocherlab_ids_gff:
    input:
        config["gff_input"],
    output:
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
    log:
        f"logs/add-kocherlab-ids-gff/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.log",
    params:
        prefix=config["id_prefix"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.6.1"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "add-prefixed-id --in-file {input} --out-file {output} --format gff --prefix {params.prefix} --attr-name kocherlab_id &> {log}"


rule add_kocherlab_ids_gtf:
    input:
        config["gtf_input"],
    output:
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gtf",
    log:
        f"logs/add-kocherlab-ids-gtf/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.log",
    params:
        prefix=config["id_prefix"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.6.1"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "add-prefixed-id --in-file {input} --out-file {output} --format gtf --prefix {params.prefix} --attr-name kocherlab_id &> {log}"


rule gff_to_sequences:
    input:
        gff_file=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        assembly_fasta=config["assembly_fasta"],
    output:
        transcripts=temp(
            f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa.tmp"
        ),
        proteins=temp(
            f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa.tmp"
        ),
    log:
        f"logs/gffread/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.log",
    singularity:
        "docker://aewebb/gffread:v0.12.7"
    resources:
        mem_mb=8000,
    threads: 1
    shell:
        "gffread -w {output.transcripts} -y {output.proteins} -g {input.assembly_fasta} {input.gff_file} &> {log}"


rule update_headers:
    input:
        gff_file=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff",
        fasta_file=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_{{seq_type}}.fa.tmp",
    output:
        f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_{{seq_type}}.fa",
    wildcard_constraints:
        seq_type="trans|pep",
    log:
        f"logs/update-header/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.{{seq_type}}.log",
    singularity:
        "docker://aewebb/pipemake_utils:v1.6.1"
    resources:
        mem_mb=2000,
    threads: 1
    shell:
        "update-header --gff-file {input.gff_file} --fasta-file {input.fasta_file} --out-file {output} --id-attr kocherlab_id &> {log}"

use rule annotations_compleasm as link_annotations_compleasm with:
    input:
        transcript_fasta=f"Processed/{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa",




