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
		"docker://huangnengcsu/compleasm:v0.2.6"
	resources:
		mem_mb=2000
	threads: 1
	shell:
		"compleasm download {params.busco_db} --library_path {params.downloads_dir}"