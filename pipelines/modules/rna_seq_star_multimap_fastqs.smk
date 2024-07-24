ruleorder: star_pair_end > star_single_end

rule all:
	input:
		expand(os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.out.bam"), sample=config['samples'])

rule star_genome_generate_rnaseq:
	input:
		fasta_file=os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa"),
		gtf_file=os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.gtf")
	output:
		index_file=os.path.join(config['paths']['index_dir'], "STAR", "SAindex")
	params:
		index_dir=directory(os.path.join(config['paths']['index_dir'], "STAR")),
		read_len=config['read_len']
	singularity:
		"docker://aewebb/star:v2.7.11b"
	resources:
		mem_mb=3200
	threads: 12
	shell:
		"""
		let "index_mem_b={resources.mem_mb} * 10**6"
		let "index_read_len={params.read_len} - 1"
		STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir} --genomeFastaFiles {input.fasta_file} --limitGenomeGenerateRAM $index_mem_b --genomeSAindexNbases 13 --sjdbGTFfile {input.gtf_file} --sjdbOverhang $index_read_len
		"""

rule star_single_end:
	input:
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		index_file=os.path.join(config['paths']['index_dir'], "STAR", "SAindex")
	output:
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Unmapped.out.mate1.gz"),
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Log.final.out")
	params:
		index_dir=os.path.join(config['paths']['index_dir'], "STAR"),
		bam_prefix=os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.")
	singularity:
		"docker://aewebb/star:v2.7.11b"
	resources:
		mem_mb=16000
	threads: 4
	shell:
		"""
		STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outReadsUnmapped Fastx --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 1 --outSAMmode None --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads}
		gzip {params.bam_prefix}Unmapped.out.mate1
		"""

rule star_pair_end:
	input:
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		r2_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R2.fq.gz"),
		index_file=os.path.join(config['paths']['index_dir'], "STAR", "SAindex")
	output:
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Unmapped.out.mate1.gz"),
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Unmapped.out.mate2.gz"),
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Log.final.out")
	params:
		index_dir=os.path.join(config['paths']['index_dir'], "STAR"),
		bam_prefix=os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.")
	singularity:
		"docker://aewebb/star:v2.7.11b"
	resources:
		mem_mb=16000
	threads: 4
	shell:
		"""
		STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outReadsUnmapped Fastx --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads}
		gzip {params.bam_prefix}Unmapped.out.mate1
		gzip {params.bam_prefix}Unmapped.out.mate2
		"""
