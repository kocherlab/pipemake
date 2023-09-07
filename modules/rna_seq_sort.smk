rule config:
	params:
		samples
		paths: 
			rnaseq_aligned_bam_dir
			rnaseq_sorted_bam_dir

rule all:
	input:
		expand(os.path.join(config['paths']['rnaseq_sorted_bam_dir'], "{sample}.sortedByCoord.bam"), sample=config['samples'])

rule sort_bam_rnaseq:
	input:
		os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.bam")
	output:
		os.path.join(config['paths']['rnaseq_sorted_bam_dir'], "{sample}.sortedByCoord.bam")
	singularity:
		"docker://biocontainers/samtools:v1.3_cv3"
	threads: 4
	shell:
		"samtools sort -@ {threads} -o {output} {input}"