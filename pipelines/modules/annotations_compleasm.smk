module config:
	params:
		species
		assembly_version
		annotation_version
		busco_database
	paths:
		annotations_dir
		downloads_dir

rule all:
	input:
		os.path.join(config['paths']['annotations_dir'], 'compleasm', 'summary.txt')

rule annotations_compleasm:
	input:
		transcript_fasta=os.path.join(config['paths']['annotations_dir'], f"{config['species']}_OGS_{config['assembly_version']}.{config['annotation_version']}_trans.fa"),
		compleasm_library=os.path.join(config['paths']['downloads_dir'], 'compleasm', f"{config['busco_database']}.done")
	output:
		os.path.join(config['paths']['annotations_dir'], 'compleasm', 'summary.txt')
	params:
		output_dir=os.path.join(config['paths']['annotations_dir'], 'compleasm'),
		download_dir=os.path.join(config['paths']['downloads_dir'], 'compleasm'),
		busco_db=config['busco_database']
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/compleasm.sif"
	threads: 20
	shell:
		"compleasm run -t{threads} -L {params.download_dir} -l {params.busco_db} -a {input.transcript_fasta} -o {params.output_dir}"
