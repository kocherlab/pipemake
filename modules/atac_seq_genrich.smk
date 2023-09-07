rule config:
	params:
		samples
		paths:
			atac_seq_sorted_bam_dir
			atac_seq_peaks_dir

rule all:
	input:
		expand(os.path.join(config['paths']['atac_seq_sorted_bam_dir'], "{sample}.sortedByCoord.bam"), sample=config['samples'])

rule sort_bam:
	input:
		os.path.join(config['paths']['atac_seq_sorted_bam_dir'], "{sample}.sortedByCoord.bam")
	output:
		os.path.join(config['paths']['atac_seq_peaks_dir'], "{sample}.narrowPeak")
	singularity:
		"docker://biocontainers/samtools:v1.3_cv3"
	threads: 4
	shell:
		"Genrich -t {input} -o {output} -v"