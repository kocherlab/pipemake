rule all:
    input:
        snail_png=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["figures_dir"],
            f"{config['species']}_{config['assembly_version']}_scaff_snail.png",
        ),
        blob_png=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["figures_dir"],
            f"{config['species']}_{config['assembly_version']}_scaff_blob.png",
        ),
        cumulative_png=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["figures_dir"],
            f"{config['species']}_{config['assembly_version']}_scaff_cumulative.png",
        ),
        os.path.join(
            config["paths"]["workflow_prefix"],
            f"{config['paths']['blobtools_dir']}_blobblurbout.tsv",
        ),

rule hifi_align_minimap2:
    input:
        hifi_fastq=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hifi_fastq_dir"],
            f"{config['hifi_sample']}.fq.gz",
        ),
        assembly_fasta=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hifi_bam_dir"],
            f"{config['hifi_sample']}.reads.bam",
        ),
    singularity:
        "docker://aewebb/minimap2:v2.28"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        "minimap2 -ax map-hifi -t {threads} {input.assembly_fasta} {input.hifi_fastq} | samtools sort -@{threads} --threads {threads} -bh -o {output}"


rule blastn_assembly_nt:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blast_dir"],
            'blastn',
            f"{config['species']}_{config['assembly_version']}.out",
        ),
    params:
        ncbi_nt_db=config["ncbi_nt_db"],
    singularity:
        "docker://aewebb/ncbi-blast:v2.16.0"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        'blastn -query {input} -db {params.ncbi_nt_db} -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads {threads} -out {output}'


rule short_assembly_records_for_blastx:
    input: 
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["blast_dir"],
                'blastx',
                f"{config['species']}_{config['assembly_version']}.short_records.fa",
            )
        ),
    params:
        max_length=config["blastx_max_length"],
    singularity:
        "docker://aewebb/seqkit:v2.10.0"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        """
        let "max_length={params.max_length} - 1"
        seqkit seq -M $max_length {input} > {output}
        """

rule blastx_short_assembly_records_diamond:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blast_dir"],
            'blastx',
            f"{config['species']}_{config['assembly_version']}.short_records.fa",
        ),
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blast_dir"],
            'blastx',
            f"{config['species']}_{config['assembly_version']}.out",
        ),
    params:
        uniprot_db=config["uniprot_db"],
    singularity:
        "docker://aewebb/diamond:v2.1.11"
    resources:
        mem_mb=56000,
    threads: 16
    shell:
        'diamond blastx --query {input} --db {params.uniprot_db} --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads {threads} > {output}'

rule blobtk_blobtools_create:
    input:
        os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.fa",
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["blobtools_dir"],
                '.create.chk'
            )
        ),
    params:
        blob_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
        )
    singularity:
        "docker://genomehubs/blobtoolkit:4.4.5"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobtools create --fasta {input} {params.blob_dir} && touch {output}"

rule blobtk_blobtools_add_cov:
    input:
        hifi_bam=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["hifi_bam_dir"],
            f"{config['hifi_sample']}.reads.bam",
        ),
        create_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
            '.create.chk'
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["blobtools_dir"],
                '.cov.chk'
            )
        ),
    params:
        blob_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
        ),
    singularity:
        "docker://genomehubs/blobtoolkit:4.4.5"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobtools add --cov {input.hifi_bam} {params.blob_dir} && touch {output}"

rule blobtk_blobtools_add_hits:
    input:
        blastn_hits=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blast_dir"],
            'blastn',
            f"{config['species']}_{config['assembly_version']}.out",
        ),
        blastx_hits=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blast_dir"],
            'blastx',
            f"{config['species']}_{config['assembly_version']}.out",
        ),
        cov_check=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
            '.cov.chk'
        )
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["blobtools_dir"],
                '.hits.chk'
            )
        ),
    params:
        blob_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
        ),
        ncbi_taxa_db=config["ncbi_taxa_db"],
    singularity:
        "docker://genomehubs/blobtoolkit:4.4.5"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobtools add --hits {input.blastn_hits} --hits {input.blastx_hits} --taxrule bestsumorder --taxdump {params.ncbi_taxa_db} {params.blob_dir} && touch {output}"

rule blobtk_blobtools_add_busco:
    input:
        busco_table=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.busco.full_table.tsv",
        ),
        cov_check=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
            '.cov.chk'
        ),
    output:
        temp(
            os.path.join(
                config["paths"]["workflow_prefix"],
                config["paths"]["blobtools_dir"],
                '.busco.chk'
            )
        ),${assembly}_blobblurbout.tsv
    singularity:
        "docker://genomehubs/blobtoolkit:4.4.5"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobtools add --busco {input.busco_table} {params.blob_dir} && touch {output}"

rule blobblurb:
    input:
        busco_table=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["assembly_dir"],
            f"{config['species']}_{config['assembly_version']}.busco.full_table.tsv",
        ),
        hits_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
            '.hits.chk'
        ),
        busco_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
            '.busco.chk'
        )
    output:
        os.path.join(
            config["paths"]["workflow_prefix"],
            f"{config['paths']['blobtools_dir']}_blobblurbout.tsv",
        )
    params:
        blob_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
        ),
    singularity:
        "docker://aewebb/blobblurb:05152024"
    resources:
        mem_mb=16000,
    threads: 1
    shell:
        "blobblurb {params.blob_dir} {input.busco_table}"

rule blobbtk_plot:
    input:
        hits_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
            '.hits.chk'
        ),
        busco_chk=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
            '.busco.chk'
        )
    output:
        snail_png=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["figures_dir"],
            f"{config['species']}_{config['assembly_version']}_scaff_snail.png",
        ),
        blob_png=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["figures_dir"],
            f"{config['species']}_{config['assembly_version']}_scaff_blob.png",
        ),
        cumulative_png=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["figures_dir"],
            f"{config['species']}_{config['assembly_version']}_scaff_cumulative.png",
        ),
    params:
        blob_dir=os.path.join(
            config["paths"]["workflow_prefix"],
            config["paths"]["blobtools_dir"],
        ),
    singularity:
        "docker://genomehubs/blobtoolkit:4.4.5"
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