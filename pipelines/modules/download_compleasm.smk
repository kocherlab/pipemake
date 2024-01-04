module config:
	params:
		species
		assembly_version
		busco_database
	paths: 
		assembly_dir
		downloads_dir

rule all:
	input:
		os.path.join(config['paths']['downloads_dir'], 'compleasm', f"{config['busco_database']}.done")

rule download_compleasm_library:
	input:
		os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa")
	output:
		os.path.join(config['paths']['downloads_dir'], 'compleasm', f"{config['busco_database']}.done")
	params:
		busco_db=config['busco_database'],
		downloads_dir=os.path.join(config['paths']['downloads_dir'], 'compleasm')
	singularity:
		"/Genomics/argo/users/aewebb/.local/images/compleasm.sif"
	shell:
		"compleasm download {params.busco_db} --library_path {params.downloads_dir}"