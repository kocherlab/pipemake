ruleorder: hisat2_pair_end > hisat2_single_end

rule all:
	input:
		expand(os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.bam"), sample=config['samples'])

rule hisat2_build:
	input:
		os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa")
	output:
		os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index.1.ht2'),
		os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index.2.ht2'),
		os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index.3.ht2'),
		os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index.4.ht2'),
		os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index.5.ht2'),
		os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index.6.ht2'),
		os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index.7.ht2'),
		os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index.8.ht2')
	params:
		index_prefix=os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index')
	singularity:
		"/Genomics/kocherlab/lab/Pipelines/imageskocherSEQ.sif"
	resources:
		mem_mb=16000
	threads: 4
	shell:
		"hisat2-build -p {threads} {input} {params.index_prefix}"

rule hisat2_single_end:
	input:
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		index1=os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index.1.ht2')
	output:
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.bam")
	params:
		index_prefix=os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index')
	singularity:
		"/Genomics/kocherlab/lab/Pipelines/imageskocherSEQ.sif"
	resources:
		mem_mb=16000
	threads: 4
	shell:
		"hisat2 --threads {threads} --dta -q -x {params.index_prefix} -U {input.r1_reads} | samtools view -@ {threads} -bh -o {output}"

rule hisat2_pair_end:
	input:
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		r2_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R2.fq.gz"),
		index1=os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index.1.ht2')
	output:
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.bam")
	params:
		index_prefix=os.path.join(config['paths']['index_dir'], 'hisat2', 'hisat2_index')
	singularity:
		"/Genomics/kocherlab/lab/Pipelines/imageskocherSEQ.sif"
	resources:
		mem_mb=16000
	threads: 4
	shell:
		"hisat2 --threads {threads} --dta -q -x {params.index_prefix} -1 {input.r1_reads} -2 {input.r2_reads} | samtools view -@ {threads} -bh -o {output}"