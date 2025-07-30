rule all:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["figures_dir"],
            f"{config['species']}_{config['assembly_version']}_scaff_snail.png",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["figures_dir"],
            f"{config['species']}_{config['assembly_version']}_scaff_blob.png",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["figures_dir"],
            f"{config['species']}_{config['assembly_version']}_scaff_cumulative.png",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            f"{config['paths']['blobtools_dir']}_blobblurbout.tsv",
        ),
        expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["blast_dir"],
                "{blast_type}",
                f"{config['species']}_{config['assembly_version']}.out",
            ),
            blast_type=["blastn", "blastx"],
        ),
#        f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_snail.png",
#        f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_blob.png",
#        f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_cumulative.png",
#        f"Assembly/blobtools/{config['species']}_{config['assembly_version']}_blobblurbout.tsv",


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
        "minimap2 -ax map-hifi -t {threads} {input.assembly_fasta} {input.hifi_fastq} | samtools sort -@{threads} --threads {threads} -O bam -o {output}"


checkpoint split_assembly:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        directory(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["assembly_dir"],
                "split",
            )
        ),
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.7"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "split-fasta --input-fasta {input} --output-dir {output}"

rule chunk_file:
    input:
        os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["assembly_dir"],
                "split",
                "{chrom}.fasta",
            )
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["assembly_dir"],
                "chunked",
                "{chrom}.chunked.fasta",
            )
        ),
#        f"BLAST/Assembly/blastn/{config['species']}_{config['assembly_version']}.out",
    params:
        chunk_size=config["chunk_size"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.2"
    resources:
        mem_mb=4000,
    threads: 1
    shell:
        "chunk-fasta --input-fasta {input} --chunk-size {params.chunk_size} --output-fasta {output}"


rule blastn_chunked_assembly_nt:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            "chunked",
            "{chrom}.chunked.fasta",
        )
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["blast_dir"],
                "blastn",
                "chunked",
                "{chrom}.chunked.out",
            )
        )
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
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            "chunked",
            "{chunk}.chunked.fasta",
        )
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["blast_dir"],
                "blastx",
                "chunked",
                "{chunk}.chunked.out",
            )
        )
#        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
#    output:
#        temp(
#            f"BLAST/Assembly/blastx/{config['species']}_{config['assembly_version']}.short_records.fa"
#        ),
    params:
        uniprot_db=config["uniprot_db"],
    singularity:
        "docker://aewebb/diamond:v2.1.11"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        "diamond blastx --query {input} --db {params.uniprot_db} --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads {threads} > {output}"

def aggregate_blast (wildcards):
    checkpoint_output = checkpoints.split_assembly.get(
        **wildcards
    ).output[0]
    return {
        "blastn": expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["blast_dir"],
                "blastn",
                "chunked",
                "{chrom}.chunked.out",
            ),
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.fasta",
                )
            ).chrom,
        ),
        "blastx": expand(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["blast_dir"],
                "blastx",
                "chunked",
                "{chrom}.chunked.out",
            ),
            chrom=glob_wildcards(
                os.path.join(
                    checkpoint_output,
                    "{chrom}.fasta",
                )
            ).chrom,
        ),
    }

rule cat_blast:
    input:
        unpack(aggregate_blast),
    output:
        blastn=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blast_dir"],
            "blastn",
            f"{config['species']}_{config['assembly_version']}.chunked.out",
        ),
        blastx=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blast_dir"],
            "blastx",
            f"{config['species']}_{config['assembly_version']}.chunked.out",
        ),
    threads: 1
    shell:
        """
        cat {input.blastn} > {output.blastn}
        cat {input.blastx} > {output.blastx}
        """

rule unchunk_blast:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blast_dir"],
            "{blast_type}",
            f"{config['species']}_{config['assembly_version']}.chunked.out",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blast_dir"],
            "{blast_type}",
            f"{config['species']}_{config['assembly_version']}.out",
        ),
=======
#        f"BLAST/Assembly/blastx/{config['species']}_{config['assembly_version']}.short_records.fa",
#    output:
#        f"BLAST/Assembly/blastx/{config['species']}_{config['assembly_version']}.out",
#    params:
#        uniprot_db=config["uniprot_db"],
    singularity:
        "docker://aewebb/pipemake_utils:v1.2.2"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "unchunk-blast --chunked-blast {input} --unchunked-blast {output}"


rule blobtk_blobtools_create:
    input:
        f"Assembly/{config['species']}_{config['assembly_version']}.fa",
    output:
        temp("Assembly/blobtools/.create.chk"),
    params:
        blob_dir=f"Assembly/blobtools/{config['species']}_{config['assembly_version']}",
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
        create_chk="Assembly/blobtools/.create.chk",
    output:
        temp("Assembly/blobtools/.cov.chk"),
    params:
        blob_dir=f"Assembly/blobtools/{config['species']}_{config['assembly_version']}",
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
        cov_check="Assembly/blobtools/.cov.chk",
    output:
        temp("Assembly/blobtools/.hits.chk"),
    params:
        blob_dir=f"Assembly/blobtools/{config['species']}_{config['assembly_version']}",
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
        cov_check="Assembly/blobtools/.cov.chk",
    output:
        temp("Assembly/blobtools/.busco.chk"),
    params:
        blob_dir=f"Assembly/blobtools/{config['species']}_{config['assembly_version']}",
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
        hits_chk="Assembly/blobtools/.hits.chk",
        busco_chk="Assembly/blobtools/.busco.chk",
    output:
        f"Assembly/blobtools/{config['species']}_{config['assembly_version']}_blobblurbout.tsv",
    params:
        blob_dir=f"Assembly/blobtools/{config['species']}_{config['assembly_version']}",
    singularity:
        "docker://aewebb/blobblurb:05152024"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobblurb {params.blob_dir} {input.busco_table}"


rule blobbtk_plot:
    input:
        hits_chk="Assembly/blobtools/.hits.chk",
        busco_chk="Assembly/blobtools/.busco.chk",
    output:
        snail_png=f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_snail.png",
        blob_png=f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_blob.png",
        cumulative_png=f"Figures/blobtools/{config['species']}_{config['assembly_version']}_scaff_cumulative.png",
    params:
        blob_dir=f"Assembly/blobtools/{config['species']}_{config['assembly_version']}",
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
