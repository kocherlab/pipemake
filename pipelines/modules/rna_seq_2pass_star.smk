ruleorder: star_pair_end_p1 > star_single_end_p1
ruleorder: star_pair_end_p2 > star_single_end_p2

rule all:
	input:
		expand(os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.bam"), sample=config['samples'])

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
		mem_mb=32000
	threads: 12
	shell:
		"""
		let "index_mem_b={resources.mem_mb} * 10**6"
		let "index_read_len={params.read_len} - 1"
		STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir} --genomeFastaFiles {input.fasta_file} --limitGenomeGenerateRAM $index_mem_b --genomeSAindexNbases 13 --sjdbGTFfile {input.gtf_file} --sjdbOverhang $index_read_len
		"""

rule star_single_end_p1:
	input:
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		index_file=os.path.join(config['paths']['index_dir'], "STAR", "SAindex")
	output:
		os.path.join(config['paths']['rnaseq_splice_aligned_dir'], "{sample}.SJ.out.tab"),
		os.path.join(config['paths']['rnaseq_splice_aligned_dir'], "{sample}.Log.final.out")
	params:
		index_dir=os.path.join(config['paths']['index_dir'], "STAR"),
		sj_prefix=os.path.join(config['paths']['rnaseq_splice_aligned_dir'], "{sample}.")
	singularity:
		"docker://aewebb/star:v2.7.11b"
	resources:
		mem_mb=16000
	threads: 4
	shell:
		"STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMmode None --outFileNamePrefix {params.sj_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads}"

rule star_pair_end_p1:
	input:
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		r2_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R2.fq.gz"),
		index_file=os.path.join(config['paths']['index_dir'], "STAR", "SAindex")
	output:
		os.path.join(config['paths']['rnaseq_splice_aligned_dir'], "{sample}.SJ.out.tab"),
		os.path.join(config['paths']['rnaseq_splice_aligned_dir'], "{sample}.Log.final.out")
	params:
		index_dir=os.path.join(config['paths']['index_dir'], "STAR"),
		sj_prefix=os.path.join(config['paths']['rnaseq_splice_aligned_dir'], "{sample}.")
	singularity:
		"docker://aewebb/star:v2.7.11b"
	resources:
		mem_mb=16000
	threads: 4
	shell:
		"STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --outSAMtype None --outFileNamePrefix {params.sj_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads}"

rule filter_star_sj_file:
	input:
		os.path.join(config['paths']['rnaseq_splice_aligned_dir'], "{sample}.SJ.out.tab")
	output:
		os.path.join(config['paths']['rnaseq_splice_aligned_dir'], "{sample}.SJ.filtered.tab")
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"cat {input} | awk '($5 > 0 && $7 > 2 && $6==0)' | cut -f1-6 | sort | uniq > {output}"

rule star_single_end_p2:
	input:
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		index_file=os.path.join(config['paths']['index_dir'], "STAR", "SAindex"),
		sj_file=os.path.join(config['paths']['rnaseq_splice_aligned_dir'], "{sample}.SJ.filtered.tab")
	output:
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.bam"),
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Log.final.out")
	params:
		index_dir=os.path.join(config['paths']['index_dir'], "STAR"),
		bam_prefix=os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.")
	singularity:
		"docker://aewebb/star:v2.7.11b"
	resources:
		mem_mb=16000
	threads: 4
	resources:
		mem_mb=1
	shell:
		"STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --sjdbFileChrStartEnd {input.sj_file} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} && "
		"mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam"
 
rule star_pair_end_p2:
	input:
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		r2_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R2.fq.gz"),
		index_file=os.path.join(config['paths']['index_dir'], "STAR", "SAindex"),
		sj_file=os.path.join(config['paths']['rnaseq_splice_aligned_dir'], "{sample}.SJ.filtered.tab")
	output:
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.bam"),
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
		"STAR --runThreadN {threads} --runMode alignReads --genomeDir {params.index_dir} --sjdbFileChrStartEnd {input.sj_file} --outSAMtype BAM Unsorted --outFileNamePrefix {params.bam_prefix} --readFilesCommand zcat --readFilesIn {input.r1_reads} {input.r2_reads} && "
		"mv {params.bam_prefix}Aligned.out.bam {params.bam_prefix}Aligned.bam"