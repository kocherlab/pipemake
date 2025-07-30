rule all:
    input:
        "Indices/STAR/SAindex",


rule star_genome_generate_rnaseq:
    input:
        fasta_file=f"Assembly/{config['species']}_{config['assembly_version']}.fa",
        gtf_file=f"Assembly/{config['species']}_{config['assembly_version']}.gtf",
    output:
        index_file="Indices/STAR/SAindex",
    params:
        index_dir="Indices/STAR",
        read_len=config["read_len"],
    singularity:
        "docker://aewebb/star:v2.7.11b"
    resources:
        mem_mb=32000,
    threads: 12
    shell:
        """
        let "index_mem_b={resources.mem_mb} * 10**6"
        let "index_read_len={params.read_len} - 1"
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir} --genomeFastaFiles {input.fasta_file} --limitGenomeGenerateRAM $index_mem_b --genomeSAindexNbases 13 --sjdbGTFfile {input.gtf_file} --sjdbOverhang $index_read_len
        """
