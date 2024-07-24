rule all:
	input:
		os.path.join(config['paths']['models_dir'], config['model_name'], f"{config['pop_name']}.pop")

rule pop_ind_file:
	input:
		os.path.join(config['paths']['models_dir'], f"{config['species']}.model")
	output:
		os.path.join(config['paths']['models_dir'], config['model_name'], f"{config['pop_name']}.pop")
	params:
		model_name=config['model_name'],
		out_dir=os.path.join(config['paths']['models_dir'], config['model_name'])
	resources:
		mem_mb=2000
	threads: 1
	singularity:
		"docker://aewebb/pipemake_utils:v0.1.27"
	shell:
		"model-pop-files --model-file {input} --model-name {params.model_name} --out-dir {params.out_dir}"
