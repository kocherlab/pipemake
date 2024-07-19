rule all:
	input:
		os.path.join(config['paths']['annotations_dir'], 'BUSCO', 'summary.txt')

rule annotations_compleasm:
	input:
		os.path.join(config['paths']['annotations_dir'], f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_pep.fa")
	output:
		os.path.join(config['paths']['annotations_dir'], 'BUSCO', 'summary.txt')
	params:
		download_dir=os.path.join(config['paths']['downloads_dir'], 'BUSCO'),
		output_dir=os.path.join(config['paths']['annotations_dir'], 'BUSCO'),
		busco_db=config['busco_database']
	singularity:
		"/Genomics/kocherlab/lab/Pipelines/images/busco.sif"
	resources:
		mem_mb=12000
	threads: 12
	shell:
		r"""
		busco -i {input} -o {params.output_dir} -l {params.busco_db} -m proteins -c {threads} --download_path {params.download_dir} -f &&
		summary_prefix=$(ls {params.output_dir} | grep -o 'short_summary.*BUSCO' | uniq) &&
		cp {params.output_dir}/${{summary_prefix}}.txt {params.output_dir}/summary.txt
		"""

