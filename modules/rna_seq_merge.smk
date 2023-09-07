rule config:
	params:
		samples
		species
		assembly_version
		paths:
			rnaseq_sorted_bam_dir
			rnaseq_bam_dir

rule all:
	input:
		os.path.join(config['paths']['rnaseq_bam_dir'], f"{config['species']}_{config['assembly_version']}.bam")

rule merge_bam:
	input:
		expand(os.path.join(config['paths']['rnaseq_sorted_bam_dir'], "{sample}.sortedByCoord.bam"), sample=config['samples'])
	output:
		os.path.join(config['paths']['rnaseq_bam_dir'], f"{config['species']}_{config['assembly_version']}.bam")
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/kocherSEQ.sif"
	threads: 4
	shell:
		"samtools merge -@ {threads} -r {output} {input}"