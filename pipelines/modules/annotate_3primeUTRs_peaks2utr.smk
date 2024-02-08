rule all:
	input:
		os.path.join(config['paths']['annotations_dir'], f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gff3"),
		os.path.join(config['paths']['annotations_dir'], f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gtf")

rule annotate_3prime_utrs_peaks2utr:
	input:
		merged_bam=os.path.join(config['paths']['rnaseq_bam_dir'], f"{config['species']}_{config['assembly_version']}.bam"),
		base_gff=os.path.join(config['paths']['annotations_dir'], f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.gff3")
	output:
		os.path.join(config['paths']['annotations_dir'], 'peaks2utr', f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gff3"),
		os.path.join(config['paths']['annotations_dir'], 'peaks2utr', f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gtf")
	params:
		out_prefix=os.path.join(config['paths']['annotations_dir'], 'peaks2utr', f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}")
	resources:
		mem_mb=120000
	threads: 10
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/peaks2utr.sif"
	shell:
		"peaks2utr {input.merged_bam} {input.base_gff} -o {params.out_prefix} -p {threads}"

rule process_peaks2utr_utrs:
	input:
		gff=os.path.join(config['paths']['annotations_dir'], 'peaks2utr', f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gff3"),
		gtf=os.path.join(config['paths']['annotations_dir'], 'peaks2utr', f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gtf")
	output:
		gff=os.path.join(config['paths']['annotations_dir'], f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gff3"),
		gtf=os.path.join(config['paths']['annotations_dir'], f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}.{config['utr_version']}.gtf")
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"""
		cp {input.gtf} {output.gtf}
		grep -v "##sequence-region" {input.gff} | grep -v "###" > {output.gff}
		"""
