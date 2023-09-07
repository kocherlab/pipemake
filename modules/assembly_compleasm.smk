rule config:
	params:
		species
		assembly_version
		busco_database
		paths:
			assembly_dir
			downloads_dir

rule all:
	input:
		os.path.join(config['paths']['assembly_dir'], 'compleasm', 'summary.txt')

rule assembly_compleasm:
	input:
		assembly_fasta=os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa"),
		compleasm_library=os.path.join(config['paths']['downloads_dir'], 'compleasm', f"{config['busco_database']}.done")
	output:
		os.path.join(config['paths']['assembly_dir'], 'compleasm', 'summary.txt')
	params:
		output_dir=os.path.join(config['paths']['assembly_dir'], 'compleasm'),
		download_dir=os.path.join(config['paths']['downloads_dir'], 'compleasm'),
		busco_db=config['busco_database']
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/compleasm.sif"
	threads: 20
	shell:
		"compleasm run -t{threads} -L {params.download_dir} -l {params.busco_db} -a {input.assembly_fasta} -o {params.output_dir}"
