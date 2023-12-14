module config:
	params:
		samples
		species
		assembly_version
	paths:
		assembly_dir
		index_dir
		rnaseq_fastq_dir
		rnaseq_aligned_bam_dir

ruleorder: star_pair_end > star_single_end

rule all:
	input:
		expand(os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.bam"), sample=config['samples'])

rule star_genome_generate_rnaseq:
	input:
		os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa")
	output:
		index_file=os.path.join(config['paths']['index_dir'], "STAR", "SAindex")
	params:
		index_dir=directory(os.path.join(config['paths']['index_dir'], "STAR"))
	resources:
		mem_mb=32000
	singularity:
		"docker://quay.io/biocontainers/star:2.7.8a--0"
	threads: 4
	shell:
		"""
		let "index_mem_b={resources.mem_mb} * 10**6"
		STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir} --genomeFastaFiles {input} --limitGenomeGenerateRAM $index_mem_b --genomeSAindexNbases 13
		"""

rule star_single_end_rnaseq:
	input:
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		index_file=os.path.join(config['paths']['index_dir'], "STAR", "SAindex")
	output:
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.bam"),
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Log.final.out")
	params:
		index_dir=os.path.join(config['paths']['index_dir'], "STAR"),
		bam_prefix=os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.")
	singularity:
		"docker://quay.io/biocontainers/star:2.7.8a--0"
	threads: 4
	shell:
		"STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} && "
		"mv {params.bam_prefix}.Aligned.out.bam {params.bam_prefix}.Aligned.bam"
 
rule star_pair_end_rnaseq:
	input:
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		r2_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R2.fq.gz"),
		index_file=os.path.join(config['paths']['index_dir'], "STAR", "SAindex")
	output:
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.bam"),
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Log.final.out")
	params:
		index_dir=os.path.join(config['paths']['index_dir'], "STAR"),
		bam_prefix=os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.")
	singularity:
		"docker://quay.io/biocontainers/star:2.7.8a--0"
	threads: 4
	shell:
		"STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads} && "
		"mv {params.bam_prefix}.Aligned.out.bam {params.bam_prefix}.Aligned.bam"