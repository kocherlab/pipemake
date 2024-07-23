rule fastp_pair_end:
	input:
		r1_reads=os.path.join(config['paths']['unfiltered_fastq_dir'], "{sample}_R1.fq.gz"),
		r2_reads=os.path.join(config['paths']['unfiltered_fastq_dir'], "{sample}_R2.fq.gz")
	output:
		r1_reads=os.path.join(config['paths']['filtered_fastq_dir'], "{sample}_R1.fq.gz"),
		r2_reads=os.path.join(config['paths']['filtered_fastq_dir'], "{sample}_R2.fq.gz"),
		json=os.path.join(config['paths']['filtered_fastq_dir'], "{sample}.json"),
		html=os.path.join(config['paths']['filtered_fastq_dir'], "{sample}.html")
	params:
		sample_prefix=os.path.join(config['paths']['filtered_fastq_dir'], "{sample}"),
		min_length=config['min_length']
	singularity:
		"docker://aewebb/pipemake_utils:v0.1.27"
	resources:
		mem_mb=16000
	threads: 4
	shell:
		"fastp -i {input.r1_reads} -I {input.r2_reads} -o {output.r1_reads} -O {output.r2_reads} --json {params.sample_prefix}.json --html {params.sample_prefix}.html --thread {threads} --detect_adapter_for_pe --cut_front --cut_tail --length_required {params.min_length}"