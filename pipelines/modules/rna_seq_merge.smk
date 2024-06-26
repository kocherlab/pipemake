rule all:
	input:
		os.path.join(config['paths']['rnaseq_bam_dir'], f"{config['species']}_{config['assembly_version']}.bam")

rule merge_bam:
	input:
		expand(os.path.join(config['paths']['rnaseq_sorted_bam_dir'], "{sample}.sortedByCoord.bam"), sample=config['samples'])
	output:
		os.path.join(config['paths']['rnaseq_bam_dir'], f"{config['species']}_{config['assembly_version']}.bam")
	singularity:
		"/Genomics/kocherlab/lab/Pipelines/imageskocherSEQ.sif"
	resources:
		mem_mb=16000
	threads: 4
	shell:
		"samtools merge -@ {threads} -r {output} {input}"