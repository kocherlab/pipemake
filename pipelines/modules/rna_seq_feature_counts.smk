rule all:
	input:
		os.path.join(config['paths']['rnaseq_count_dir'], f"{config['species']}_{config['assembly_version']}.featurecounts.tsv")

ruleorder: feature_counts_pair_end > feature_counts_single_end

rule feature_counts_pair_end:
	input:
		gtf_file=os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.gtf"),
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		r2_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R2.fq.gz"),
		sorted_bam=os.path.join(config['paths']['rnaseq_sorted_bam_dir'], "{sample}.sortedByCoord.bam")
	output:
		os.path.join(config['paths']['rnaseq_count_dir'], "featureCounts", "{sample}.featurecounts.txt")
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/kocherSEQ.sif"
	resources:
		mem_mb=4000
	threads: 4
	shell:
		"featureCounts -p -T {threads} -a {input.gtf_file} -o {output} {input.sorted_bam}"

rule feature_counts_single_end:
	input:
		gtf_file=os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.gtf"),
		r1_reads=os.path.join(config['paths']['rnaseq_fastq_dir'], "{sample}_R1.fq.gz"),
		sorted_bam=os.path.join(config['paths']['rnaseq_sorted_bam_dir'], "{sample}.sortedByCoord.bam")
	output:
		os.path.join(config['paths']['rnaseq_count_dir'], "featureCounts", "{sample}.featurecounts.txt")
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/kocherSEQ.sif"
	resources:
		mem_mb=4000
	threads: 4
	shell:
		"featureCounts -T {threads} -a {input.gtf_file} -o {output} {input.sorted_bam}"

rule feature_counts_report:
	input:
		expand(os.path.join(config['paths']['rnaseq_count_dir'], "featureCounts", "{sample}.featurecounts.txt"), sample=config['samples'])
	output:
		os.path.join(config['paths']['rnaseq_count_dir'], f"{config['species']}_{config['assembly_version']}.featurecounts.tsv")
	params:
		count_dir=os.path.join(config['paths']['rnaseq_count_dir'], "featureCounts")
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/BioPython.sif"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"featureCounts_report.py {params.count_dir} {output}"
