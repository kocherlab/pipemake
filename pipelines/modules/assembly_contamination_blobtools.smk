rule all:
    input:
        f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_snail.png",
        f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_blob.png",
        f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_cumulative.png",
        f"Tables/blobtools/{config['species']}_{config['assembly_version']}_blobblurbout.tsv",
        expand(
            f"BLAST/Assembly/{{blast_type}}/{config['species']}_{config['assembly_version']}.out",
            blast_type=["blastn", "blastx"],
        ),


rule hifi_align_minimap2:
    input:
        hifi_fastq=f"HiFi/FASTQ/{config['species']}.fq.gz",
        assembly_fasta=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        f"HiFi/BAM/Aligned/{config['species']}.reads.bam",
    singularity:
        "docker://aewebb/minimap2:v2.28"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        "minimap2 -ax map-hifi -t {threads} {input.assembly_fasta} {input.hifi_fastq} | samtools sort --threads {threads} -O bam -o {output}"


checkpoint split_assembly:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        directory("Assembly/split"),
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "split-fasta --input-fasta {input} --output-dir {output}"


rule chunk_file:
    input:
        "Assembly/split/{chrom}.fasta",
    output:
        temp("Assembly/chunked/{chrom}.chunked.fasta"),
    params:
        chunk_size=config["chunk_size"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "chunk-fasta --input-fasta {input} --chunk-size {params.chunk_size} --output-fasta {output}"


rule blastn_chunked_assembly_nt:
    input:
        "Assembly/chunked/{chrom}.chunked.fasta",
    output:
        temp("BLAST/Assembly/blastn/chunked/{chrom}.chunked.out"),
    params:
        ncbi_nt_db=config["ncbi_nt_db"],
    singularity:
        "docker://ncbi/blast:2.16.0"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        'blastn -query {input} -db {params.ncbi_nt_db} -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads {threads} -out {output}'


rule blastx_chunked_assembly_records_diamond:
    input:
        "Assembly/chunked/{chunk}.chunked.fasta",
    output:
        temp("BLAST/Assembly/blastx/chunked/{chunk}.chunked.out"),
    params:
        uniprot_db=config["uniprot_db"],
    singularity:
        "docker://aewebb/diamond:v2.1.11"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        "diamond blastx --query {input} --db {params.uniprot_db} --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads {threads} > {output}"


def aggregate_blast(wildcards):
    checkpoint_output = checkpoints.split_assembly.get(**wildcards).output[0]
    return {
        "blastn": expand(
            "BLAST/Assembly/blastn/chunked/{chrom}.chunked.out",
            chrom=glob_wildcards(f"{checkpoint_output}/{{chrom}}.fasta").chrom,
        ),
        "blastx": expand(
            "BLAST/Assembly/blastx/chunked/{chrom}.chunked.out",
            chrom=glob_wildcards(f"{checkpoint_output}/{{chrom}}.fasta").chrom,
        ),
    }


rule cat_blast:
    input:
        unpack(aggregate_blast),
    output:
        blastn=f"BLAST/Assembly/blastn/{config['species']}_{config['assembly_version']}.chunked.out",
        blastx=f"BLAST/Assembly/blastx/{config['species']}_{config['assembly_version']}.chunked.out",
    threads: 1
    shell:
        """
        cat {input.blastn} > {output.blastn}
        cat {input.blastx} > {output.blastx}
        """


rule unchunk_blast:
    input:
        f"BLAST/Assembly/{{blast_type}}/{config['species']}_{config['assembly_version']}.chunked.out",
    output:
        f"BLAST/Assembly/{{blast_type}}/{config['species']}_{config['assembly_version']}.out",
    singularity:
        "docker://aewebb/pipemake_utils:v1.3.2"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "unchunk-blast --chunked-blast {input} --unchunked-blast {output}"


rule blobtk_blobtools_create:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        temp("Tables/blobtools/.create.chk"),
    params:
        blob_dir=f"Tables/blobtools/{config['species']}_{config['assembly_version']}",
    singularity:
        "docker://genomehubs/blobtoolkit:4.4.6"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobtools create --fasta {input} {params.blob_dir} && touch {output}"


rule blobtk_blobtools_add_cov:
    input:
        hifi_bam=f"HiFi/BAM/Aligned/{config['species']}.reads.bam",
        create_chk="Tables/blobtools/.create.chk",
    output:
        temp("Tables/blobtools/.cov.chk"),
    params:
        blob_dir=f"Tables/blobtools/{config['species']}_{config['assembly_version']}",
    singularity:
        "docker://genomehubs/blobtoolkit:4.4.6"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobtools add --cov {input.hifi_bam} {params.blob_dir} && touch {output}"


rule blobtk_blobtools_add_hits:
    input:
        blastn_hits=f"BLAST/Assembly/blastn/{config['species']}_{config['assembly_version']}.out",
        blastx_hits=f"BLAST/Assembly/blastx/{config['species']}_{config['assembly_version']}.out",
        cov_check="Tables/blobtools/.cov.chk",
    output:
        temp("Tables/blobtools/.hits.chk"),
    params:
        blob_dir=f"Tables/blobtools/{config['species']}_{config['assembly_version']}",
        ncbi_taxa_db=config["ncbi_taxa_db"],
    singularity:
        "docker://genomehubs/blobtoolkit:4.4.6"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobtools add --hits {input.blastn_hits} --hits {input.blastx_hits} --taxrule bestsumorder --taxdump {params.ncbi_taxa_db} {params.blob_dir} && touch {output}"


rule blobtk_blobtools_add_busco:
    input:
        busco_table=f"Assembly/{config['species']}_{config['assembly_version']}.busco.full_table.tsv",
        cov_check="Tables/blobtools/.cov.chk",
    output:
        temp("Tables/blobtools/.busco.chk"),
    params:
        blob_dir=f"Tables/blobtools/{config['species']}_{config['assembly_version']}",
    singularity:
        "docker://genomehubs/blobtoolkit:4.4.6"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobtools add --busco {input.busco_table} {params.blob_dir} && touch {output}"


rule blobblurb:
    input:
        busco_table=f"Assembly/{config['species']}_{config['assembly_version']}.busco.full_table.tsv",
        hits_chk="Tables/blobtools/.hits.chk",
        busco_chk="Tables/blobtools/.busco.chk",
    output:
        f"Tables/blobtools/{config['species']}_{config['assembly_version']}_blobblurbout.tsv",
    params:
        blob_dir=f"Tables/blobtools/{config['species']}_{config['assembly_version']}",
    singularity:
        "docker://aewebb/blobblurb:05152024"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobblurb {params.blob_dir} {input.busco_table}"


rule blobbtk_plot:
    input:
        hits_chk="Tables/blobtools/.hits.chk",
        busco_chk="Tables/blobtools/.busco.chk",
    output:
        snail_png=f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_snail.png",
        blob_png=f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_blob.png",
        cumulative_png=f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_cumulative.png",
    params:
        blob_dir=f"Tables/blobtools/{config['species']}_{config['assembly_version']}",
    singularity:
        "docker://genomehubs/blobtoolkit:4.4.6"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        """
        sed -i '/level/d' {params.blob_dir}/meta.json
        blobtk plot -v snail -d {params.blob_dir} -o {output.snail_png}
        blobtk plot -v blob -d {params.blob_dir} -o {output.blob_png}
        blobtk plot -v cumulative -d {params.blob_dir} -o {output.cumulative_png}
        """
